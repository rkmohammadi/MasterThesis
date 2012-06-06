function   [lbl_bin_FRST, en_im, BinarySOpenByRec, edg]= PreProcess(I, minSize)
% function   [lbl_bin_FRST, en_im, BinarySOpenByRec]= PreProcess(I)
% this function process an image by first enhancing the image guality using
% bilateral filter. Then apply fast radial symmetry transform and some
% morphological operation to remove noisy regions. Finally label the image
% input:
%      I, either filename(string) or color image matrix in RGB. default value is 
%       '../images/GATA-4/M1G40033.TIF'
%      minSize, an integer that is minimum size of acceptable regions
%       default value is 20
% output: 
%      lbl_bin_FRST, label image in dip_image class,  
%      en_im, enhanced image after applying bilateral filter
%      BinarySOpenByRec, binary image before labeling
%      edg, edge map of the image using Canny edge detector on L-channel of
%      the enhanced image.
%
% Note: in labeling process, all regions smaller than 20 pixeles are removed.
%
% Author: A.Rahim Kadkhodamohammadi (r.k.mohammadi@gmail.com)
% 24/02/2012 CBA, Uppsala University
%--------------------------------------------------------------


% check input
switch nargin
    case 0
        I =  '../images/GATA-4/M1G40033.TIF';% default image
    case 1
        minSize = 20; % defalut minimum size of acceptable region
    case 2
        disp('Processing image ...');
    otherwise
        error('Please check function for valid inputs!!!');
end

if ischar(I)
    Im = imread(I);
else
    Im = I;
end

sz = size(Im);
% preprocessing 
doubleIm = double(Im) ./255;
w = 5;
sigma = [3, .2];
en_im = bfilter2(doubleIm, w, sigma);
llab =RGB2Lab(en_im);
edg =edge(llab(:,:,1), 'canny', [], 1 );
%%% save edge if needed
%save([I(18:end-4) 'edg.mat'], 'edg');
% fast radial symmetry of s plane in HSV 
hsv_im = rgb2hsv(en_im);

% convert s channel to uint8 
s_channel = percentileStretch(hsv_im(:,:,2), [0, 1]);

gaussK = fspecial('gaussian',[7 7],2);
s_channel = imfilter(s_channel, gaussK,'same');

% it worths to try [3 4 5 6 7]?
[S, So] = fastradial(s_channel, [3 4 5], 1, .2); %[4 5 6 7]

% figure 
% subplot(2,1,1), imshow(S), title('a) FRST of S pplane in HSV');
% since all cells are bright in S, then the positive value in FRST indication 
% nuclei/cell
bin_s = S >1; %im2bw(S,graythresh(S));
% subplot(2,1,2), imshow(bin_s), title('binary image of a');

% apply opening by reconstruction to remove small noise

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

% figure, imshow(BinarySOpenByRec), 
% title('The remaining cell after opening by reconstruction');

% to show the seed on top of the origenla image uncomment following 

% figure, imshow(Im);
% maskIm=Im; %zeros(size(Im));
% tmp = Im(:,:,2);
% tmp(BinarySOpenByRec) = 255;
% maskIm(:,:,2) = tmp; %BinarySOpenByRec.*255; %255;
% hold on
% maskH = imshow(maskIm);
% set(maskH, 'AlphaData', 0.8);
% hold off 

% label Region and colect some properties
lbl_bin_FRST = label(BinarySOpenByRec,1,minSize,0);
end