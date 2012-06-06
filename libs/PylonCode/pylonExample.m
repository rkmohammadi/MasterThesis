function pylonExample

%---------------------------------------------------------
%%setting up
vl_setup
if exist('SolveQPBO') ~= 3
    if ~exist('vl_setup')
        error('To run this example you also need to install the vl_feat toolbox.');
    end
    pylonSetup;
end

%----------------------------------------------------------
%% Creating tree hierarchy from the segmentation tree
im = imread('bsd_example.jpg');
brush = imread('bsd_example.png');
ucm = imread('bsd_ucm.bmp');

regions = bwlabel(ucm == 0,4);
edges = find(ucm > 0);
[edgeI edgeJ] = ind2sub(size(ucm),edges);
edgeStrength = ucm(edges);

regions = padarray(regions,[1 1]);
edges = sub2ind(size(regions), edgeI+1, edgeJ+1); %reflect the padding

neighbors4 = [regions(edges-1) regions(edges+1) regions(edges+size(regions,1)) regions(edges-size(regions,1))];
neighbors8 = [neighbors4 ...
    regions(edges-1-size(regions,1)) regions(edges+1-size(regions,1)) ...
    regions(edges-1+size(regions,1)) regions(edges+1+size(regions,1))];

neighbors4 = neighbors4';
neighbors8 = neighbors8';

nRegions = max(regions(:));

bLength4 = zeros(nRegions);
bLength8 = zeros(nRegions);
bStrength = zeros(nRegions);

for i = 1:size(neighbors8,2)
    t = unique(neighbors8(:,i));
    bLength8(t(2:end),t(2:end))=bLength8(t(2:end),t(2:end))+1;
    t = unique(neighbors4(:,i));
    bLength4(t(2:end),t(2:end))=bLength4(t(2:end),t(2:end))+1;   
end

bLength8(bLength8==1) = 0; %removing "accross corner" neighbors 
bLength = (bLength4+bLength8)*0.5; %reasonable approximation to the Euclidean length

bLength = bLength-diag(diag(bLength));

[nbr1 nbr2] = ind2sub(size(bLength),find(bLength));
t = find(nbr1 > nbr2);
nbr1(t) = [];
nbr2(t) = [];

where8 = cell(nRegions,1);

for i = 1:nRegions
    where8{i} = find(any(neighbors8 == i));
end

for i = 1:numel(nbr1)
    bStrength(nbr1(i),nbr2(i)) = median(edgeStrength(intersect(where8{nbr1(i)},where8{nbr2(i)})));
end
bStrength = max(bStrength,bStrength');

%filling in edge pixels
el = strel('diamond',1); 
for i = 1:2
   tmp = imdilate(regions,el);
   regions(regions == 0) = tmp(regions == 0);
end

bStrength(bStrength == 0) = +inf;
bStrength(sub2ind(size(bStrength),1:size(bStrength,1),1:size(bStrength,1))) = 0;

Tree = linkage(squareform(bStrength));
baseRegions = regions(2:end-1,2:end-1);
nBaseRegions = nRegions;

[i j s] = find(bLength);
where = i >= j;
i(where) = [];
j(where) = [];
s(where) = [];
t = bStrength(sub2ind(size(bStrength), i, j));
V = [i j s.*exp(-t/50)*10]';
V0 = [i j s.*exp(-t/50)*0]'; %disabled pairwise terms for comparison


%------------------------------------------------------------
%% computing unary terms based on color histograms
R = im(:,:,1);
G = im(:,:,2);
B = im(:,:,3);
colorData = [R(:) G(:) B(:)]';

Rb = brush(:,:,1);
Gb = brush(:,:,2);
Bb = brush(:,:,3);
brushData = [Rb(:) Gb(:) Bb(:)]';

K = 64;
vocab = vl_ikmeans(colorData, K);
words = vl_ikmeanspush(colorData, vocab);

brushes = { (Rb == 255 & Gb == 0 & Bb == 0), ...
    (Rb == 0 & Gb == 255 & Bb == 0), ...
    (Rb == 0 & Gb == 0 & Bb == 255), ...
    (Rb == 255 & Gb == 255 & Bb == 0) };


brushWords = { words(brushes{1}), words(brushes{2}), words(brushes{3}), words(brushes{4})};

brushHist = ...
    [accumarray(brushWords{1}',ones(size(brushWords{1}')),[K 1]),...
    accumarray(brushWords{2}',ones(size(brushWords{2}')),[K 1]),...
    accumarray(brushWords{3}',ones(size(brushWords{3}')),[K 1]),...
    accumarray(brushWords{4}',ones(size(brushWords{4}')),[K 1])];

%approximating chi-square kernel
brushHist = bsxfun(@times,brushHist,1.0./sum(brushHist));
brushHist = vl_homkermap(brushHist,1);

%filling in the color histograms for leaf regions
regionHist = zeros(K, nBaseRegions*2-1);
t = accumarray((baseRegions(:)-1)*K+double(words(:)),ones(numel(words(:)),1),[nBaseRegions*K 1]);
regionHist(:,1:nBaseRegions) = reshape(t, [K nBaseRegions]); 

%filling in the color histograms for non-leaf regions
for i = 1:size(Tree,1)
    regionHist(:,nBaseRegions+i) =  regionHist(:,Tree(i,1))+regionHist(:,Tree(i,2));
end
regionSizes = sum(regionHist);

%approximating chi-square kernel
regionHist = vl_homkermap(bsxfun(@times,regionHist,1./regionSizes),1);
regionHist = bsxfun(@times,regionHist,regionSizes);

%computing unary term
U = -brushHist'*regionHist;
Uflat = U;
Uflat(:,nBaseRegions+1:end) = 1e6; 

%------------------------------------------------------------
%% computing hard unaries based on brushes
hardU = zeros(4,nBaseRegions);
hardU(1,baseRegions(brushes{1})) = -1e6;
hardU(2,baseRegions(brushes{2})) = -1e6;
hardU(3,baseRegions(brushes{3})) = -1e6;
hardU(4,baseRegions(brushes{4})) = -1e6;

%------------------------------------------------------------
%% optimization
tic
xlabels_full = pylonInferenceMultiClass(nBaseRegions, Tree, U, V, hardU);
flabels_full = pylonConvertLabels(xlabels_full, Tree, nBaseRegions);
timing = toc
disp(['Time spent in the pylon inference is ' num2str(timing) ' seconds.']);

xlabels_nobnd = pylonInferenceMultiClass(nBaseRegions, Tree, U, V0, hardU);
flabels_nobnd = pylonConvertLabels(xlabels_nobnd, Tree, nBaseRegions);

xlabels_flat = pylonInferenceMultiClass(nBaseRegions, Tree, Uflat, V, hardU);
flabels_flat = pylonConvertLabels(xlabels_flat, Tree, nBaseRegions);



%------------------------------------------------------------
%% visualizing results
close all;
figure('Name','image from Berkeley segmentation dataset');
imshow(im);

figure('Name','segmentation tree (UCM method [Arbelaez et al.])');
imshow(ucm);

figure('Name','user brushes');
imshow(brush);

palette = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0; 1.0 1.0 0.0];

figure('Name','pylon result');
imshow(pylonVisualize(xlabels_full, nBaseRegions, baseRegions, Tree, palette));

figure('Name','pylon result with pairwise terms switched off');
imshow(pylonVisualize(xlabels_nobnd, nBaseRegions, baseRegions, Tree, palette));

figure('Name','flat CRF result');
imshow(pylonVisualize(xlabels_flat, nBaseRegions, baseRegions, Tree, palette));




