% manually define class of some region
% the file name are stored invariable images

images = {'../images/GATA-4/M3G4001.TIF'; ...
    '../images/GATA-4/M3G4004.TIF'; ...
    '../images/GATA-4/M3G4006.TIF'; ...
    '../images/GATA-4/M1G40019.TIF'; ...
    '../images/GATA-4/M1G40022.TIF'; ...
    '../images/GATA-4/M1G40036.TIF'};


for ii = 1 : numel(images)
    imgFile = images{ii};
    load([imgFile(18:end-4) 'WS.mat']);
    
    GMandSertoli = or(GMtuned,sertoli);
    [~, h] = overlayCluster(en_im, GMandSertoli, [], 1); 

    hold on
    % sertoliLbl = bwlabel(sertoli,4);
    stats = regionprops(sertoliLbl, 'Centroid');
    num = numel(stats);

    for i = 1 : num
        NS = closeGM{i};
        for j =1 : size((NS{2}),1)
            if isempty(NS{3}{j}) %&& NS{2}(j,5) == 0
                line([stats(i).Centroid(1), NS{2}(j,4)], ...
                    [stats(i).Centroid(2), NS{2}(j,3)]);
                break
            end
        end
    end
    
    fprintf(1,'Choose Sertoli class(class1Region) regions \n');
    [xi, yi , ~] = impixel;
    class1Region = zeros(length(xi)-1,1);
    for r=1 : length(xi)-1
        rlbl= sertoliLbl(yi(r), xi(r));
        class1Region(r) = rlbl;
    end
    
    fprintf(1, 'Choose germ class(class2Region) regions \n');
    [xi, yi , ~] = impixel;
    class2Region = zeros(length(xi)-1,1);
    for r=1 : length(xi)-1
        class2Region(r) = lbl_bin_FRST(yi(r), xi(r));
    end
    
    fprintf(1, 'Choose boarder class(class3Region) regions \n');
    [xi, yi , ~] = impixel;
    class3Region = zeros(length(xi)-1,1);
    for r=1 : length(xi)-1
        class3Region(r) = lbl_bin_FRST(yi(r), xi(r));
    end
    
    fprintf(1, 'Choose semi-Sertoli class(class4Region) regions \n');
    [xi, yi , ~] = impixel;
    class4Region = zeros(length(xi)-1,1);
    for r=1 : length(xi)-1
        class4Region(r) = lbl_bin_FRST(yi(r), xi(r));
    end
    
    save([imgFile(1:end-3) 'mat'], 'en_im', 'lbl_bin_FRST', ...
        'class1Region', 'class2Region', 'class3Region', 'class4Region');
end

%% show the selected region if you want check
for ii =1 : numel(images)
    imgFile = images{ii};
    load([imgFile(1:end-3) 'mat']);
    
    figure(92),
    imshow(en_im);
    hold on; 
    mask = zeros(size(lbl_bin_FRST));
    for j =1 : length(class2Region)
        if class2Region(j) ~= 0
            mask = mask + (lbl_bin_FRST== class2Region(j));
        end
    end
    mask = mask* 255; 
    % if we have value greater than 255 means something is wrong
    if max(mask(:)) > 255
        msgbox(['Please check the label it is impossible to ', ...
            'have region greater than 255']);
    end
    colorMask = uint8(cat(3, mask, mask, mask));
    hold on
    handel = imshow(colorMask); title(['regions on top of the image', imgFile]);
    set(handel,'AlphaData',0.4);
    hold off
end

%% apply measure to compute some properties for the region

measurmentIDs = {'Size', 'CartesianBox', 'Minimum', 'Maximum', ...
    'Perimeter', 'Feret', 'Radius', 'ConvexArea', 'ConvexPerimeter', ...
    'P2A', 'PodczeckShapes', 'Convexity', 'Sum','Mean', 'StdDev', ...
    'Skewness', 'ExcessKurtosis', 'MaxVal', 'MinVal', 'Center', 'Inertia',...
    'Mu', 'MajorAxes', 'Center'};
lab = RGB2Lab(en_im);
stretchedL = percentileStretch(lab(:,:,1), [0,1]);
diplbl = label(BinarySOpenByRec,1,20,0);
msr = measure(diplbl, mat2im(stretchedL),measurmentIDs, [],1);
% sort msr based on Id
tmp = msr.id;
[~, sortId] = sort(tmp);

%%
feildName = fieldnames(msr);
msrC1 = msr(class1Region);
msrC2 = msr(class2Region);

for fn =1 : length(feildName)
    t1 = getfield(msrC1, feildName{fn});
    t2 = getfield(msrC2, feildName{fn});
    switch size(t1,1)
        case 1
            h= scatterplot([t1; 1:length(t1)]');
            hold on
            scatterplot([t2; 1:length(t2)]',1,1,'+ red',h);
            title(feildName{fn});
            hold off
            
        case 2
            h= scatterplot(t1');
            hold on
            scatterplot(t2',1,1,'+ red',h);
            title(feildName{fn});
            hold off
        otherwise
            fprintf(1,'the size of %s is %d,%d \n',feildName{fn}, size(t1));
    end
end

%%
t1 = features(class1Region, :) * 300;
t2 = features(class2Region, :) * 300;

ph1 = [1: size(t1,1)]';
ph2 = [1: size(t2,1)]';

for f =1 : size(features,2)
            h= scatterplot([t1(:,f), ph1]);
            hold on
            scatterplot([t2(:,f), ph2],1,1,'+ red',h);
            title(num2str(f));
            hold off
end
%%
fIdx = [7,8,10];%10,20,21,22,27,28,44,45,46];
msrROI = msr; %msr([class1Region; class2Region]);
rawFeatures = [msrROI.Mean; msrROI.MinVal; msrROI.Perimeter; msrROI.Size;...
    msrROI.StdDev; msrROI.Skewness; msrROI.ExcessKurtosis]';
% rawFeatures = [msrROI.Skewness; msrROI.ExcessKurtosis]';
% rawFeatures = [rawFeatures, features(:,fIdx)];

[fr, ~]=size(rawFeatures);
nRawFeatures = (rawFeatures - ones(fr,1)*mean(rawFeatures)) ./ (ones(fr,1)*std(rawFeatures));

K = 2; %4
opts = statset('Display','final');
[idx,ctrs,sumdist] = kmeans(nRawFeatures,K,'dist','sqEuclidean', 'start', 'cluster',...
'replicates',5,'EmptyAction','singleton','Options',opts);

rId = msrROI.id;
y= zeros(size(lbl_bin_FRST));
for ii =1 : length(rId)
    y(lbl_bin_FRST == rId(ii)) = idx(ii);
end

sz = size(en_im);
for ii = 1: K
    figure,imshow(en_im);
    RGBDist = uint8(cat(3, zeros(sz(1:2)),(y == ii)*255 ,zeros(sz(1:2))));
    hold on
    distIm = imshow(RGBDist); title(['Cluster ', int2str(ii) ' out of ', int2str(K)]);
    set(distIm,'AlphaData',0.3);
end
