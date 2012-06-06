images = {'../images/GATA-4/M3G4001.TIF'; ...
    '../images/GATA-4/M3G4004.TIF'; ...
    '../images/GATA-4/M3G4006.TIF'; ...
    '../images/GATA-4/M1G40019.TIF'; ...
    '../images/GATA-4/M1G40022.TIF'; ...
    '../images/GATA-4/M1G40036.TIF'};

imgFile = images{1};
Im = imread(imgFile);
load([imgFile(1:end-3) 'mat']);

lab = RGB2Lab(Im);
for i =1 : 3
   labFull(:,:,i) = percentileStretch(lab(:,:,i), [0 1]); 
end

stats = regionprops(lbl_bin_FRST, 'Centroid', 'MajorAxisLength', ...
    'PixelList', 'BoundingBox', 'SubarrayIdx');

%% find neighbors 

% neighborhood region gragh
% this array will store the lebel of neighboring regions at the specified
% distance

neiborhood = 35;
tic
NRG = findNeighbors(lbl_bin_FRST, stats, neiborhood);
nBins = 50;
regions = RegionHist(lab, stats, nBins);
toc
%%
% THRESHOLD get a similarity of 0.
THRESHOLD= 3;
% The normalization factor. Should be 0 <= m < 1. 
% 0.9 experimentally yielded good results. 0.5 is the generalization of
% chi^2 which also yields good results.
m= 0.9;

% The two histograms P and Q

% The sparse bin-similarity matrix. See other demos for fast mex
% computation of this kind of matrix.
A= sparse(nBins,nBins);
for i=1:nBins
    for j=max([1 i-THRESHOLD+1]):min([nBins i+THRESHOLD-1])
        A(i,j)= 1-(abs(i-j)/THRESHOLD); 
    end
end

% The demo includes several ways to call QC
% demo_QC_compute(P, Q, A, m); dist1= QC_full_sparse(P,Q,A,m);
%% Check similarities of histograms
numChannel = size(regions(1).hist, 1);
for i = 1 : size(NRG)
    for ch = 1 : numChannel
        for j = 1 : size(NRG{i}, 2)
            regions(i).similarity(j,ch) = QC_full_sparse(regions(i).hist(ch,:), ...
                regions(NRG{i}(j)).hist(ch,:), A, m);
%             corr(regions(i).hist(ch,:)', ...
%                 regions(NRG{i}(j)).hist(ch,:)');
        end
    end
end

%% make feature vector

trainMatrix = [];
trainlbl = [];
k =1;
for j =1 : length(class1Region)
    if class1Region(j) ~= 0
        hists= [];
        for ch =1 : 1 % numChannel
            hists = [hists, regions(j).hist(ch, 10:49)];
        end
        trainMatrix(k,:) = hists;
        trainlbl(k) = 1;
        k = k+1;
    end
end

for j =1 : length(class2Region)
    if class2Region(j) ~= 0
        hists= [];
        for ch =1 : 1 %numChannel
            hists = [hists, regions(j).hist(ch, 10:49)];
        end
        trainMatrix(k,:) = hists;
        trainlbl(k) = -1;
        k = k+1;
    end
end

testMatrix = [];
for r = 1 : length(regions)
    hists = [];
    for ch =1 : numChannel
        hists = [hists, regions(j).hist(ch, :)];
    end
    testMatrix(r, : ) = hists;
end

model = svmtrain(trainlbl', trainMatrix, '-c .125, -g .125');
[prelbl, accuracy, prob_est] = svmpredict(ones(size(regions)), testMatrix, model);


