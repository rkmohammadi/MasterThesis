% doing some experiments by creating a groph and some cue for each region
% the try to make separate graph for each segment of grem mass based  upon 
% the cues 
% algorithm 
% - read the image and do some enhancement bilateral filter
% - compute fast radial symmetry 
% - create regions 
% - create region connectivity graph
% - cut edges that connect non grem cells 

%% read images 
addlibs; 

% '../images/GATA-4/M3G4008.TIF'  '../images/GATA-4/M1G40048.TIF'
fileName =  '../images/GATA-4/M1G40033.TIF';%M1G40049.TIF';
Im = imread(fileName);
sz = size(Im);
%% preprocessing 
doubleIm = double(Im) ./255;
w = 5;
sigma = [3, .2];
en_im = bfilter2(doubleIm, w, sigma);

%% fast radial symmetry of s plane in HSV 
hsv_im = rgb2hsv(en_im);

% convert s channel to uint8 
s_channel = percentileStretch(hsv_im(:,:,2), [0, 1]);

gaussK = fspecial('gaussian',[7 7],2);
s_channel = imfilter(s_channel, gaussK,'same');

% it worths to try with [3 4 5 6 7]?
[S, So] = fastradial(s_channel, [3 4 5], 1, .2);

% figure 
% subplot(2,1,1), imshow(S), title('a) FRST of S pplane in HSV');
% since all cells are bright in S, then the positive value in FRST indication 
% nuclei/cell
bin_s = S >1; %im2bw(S,graythresh(S));
% subplot(2,1,2), imshow(bin_s), title('binary image of a');

%% apply opening by reconstruction to remove small noise
SE3 = strel('disk', 3); % big SE to only preserve the lumen 
erodedBinaryS = imerode(bin_s, SE3);
BinarySOpenByRec = imreconstruct(erodedBinaryS, bin_s);

% dilatedBinaryS = imdilate(BinarySOpenByRec, SE3);
% BinarySOpenByRec = imreconstruct(dilatedBinaryS, BinarySOpenByRec);
SE2 = strel('disk', 2);
SE1 = strel('disk', 1);
BinarySOpenByRec = imclose(imopen(BinarySOpenByRec, SE2), SE2);

% grow to have outer boarder of the regions
BinarySOpenByRec = imdilate(BinarySOpenByRec, SE1);

figure, imshow(BinarySOpenByRec), 
title('The remaining cell after opening by reconstruction');

%% show the seed on top of the origenla image
% Blue= Im(:,:,3);
% Blue(BinarySOpenByRec) = 255;
% figure, imshow(Blue)
figure, imshow(Im);
maskIm=Im; %zeros(size(Im));
tmp = Im(:,:,2);
tmp(BinarySOpenByRec) = 255;
maskIm(:,:,2) = tmp; %BinarySOpenByRec.*255; %255;
hold on
maskH = imshow(maskIm);
set(maskH, 'AlphaData', 0.8);
hold off 
%% stretching channel
ImFull = Im;
hsvFull = hsv_im;
lab = RGB2Lab(Im);
for i =1 : 3
    ImFull(:,:,i) = percentileStretch(Im(:,:,i),[0, 1]);
    hsvFull(:,:,i) = percentileStretch(hsv_im(:,:,i),[0, 1]);
    labFull(:,:,i) = percentileStretch(lab(:,:,i), [0 1]); 
end

%% label Region and colect some properties

% [lbl_bin_FRST num] = bwlabel(BinarySOpenByRec,4);
lbl_bin_FRST = double(label(BinarySOpenByRec,1,20,0));
stats = regionprops(lbl_bin_FRST, 'Centroid', 'MajorAxisLength', ...
    'PixelList', 'BoundingBox', 'SubarrayIdx');

%% find neighbors 

% neighborhood region gragh
% this array will store the lebel of neighboring regions at the specified
% distance

neiborhood = 35;
tic
NRG = findNeighbors(lbl_bin_FRST, stats, neiborhood);
nBins = 255;
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

%% show graph
figure, imshow(en_im), %lbl_bin_FRST), 
% title(['connected regions', fileName]);
stats1 = regionprops(double(label(GCCluster == 1)),'Centroid');
for i =1 : size(stats1,1)
    NS = NRG1{i};
    for j =1 : numel(NS)
        %if (regions(i).similarity(j,1)) > 3
            line([stats1(i).Centroid(1), stats1(NS(j)).Centroid(1)], ...
                [stats1(i).Centroid(2), stats1(NS(j)).Centroid(2)]);
        %end
    end
end

%%
FM = -ones(3,3);
FM(2,2) = 8;
T = zeros(sz(1:2));
for ii=1 : 3
    T = T + abs(imfilter(lab(:,:,i), FM, 'same'));
end