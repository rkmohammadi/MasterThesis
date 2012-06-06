% try different channel to find a good color space which celles are 
% better separable 
addlibs; % add libraries 
% some images M1G40025.TIF, M1G40027.TIF
% another set M3G4001.TIF, M3G4008.TIF, M3G4009.TIF, M3G4012.TIF  
% fileName = '../images/set_1/M1x10ex1.tif';
fileName =  '../images/GATA-4/M1G40024.TIF';
Im = imread(fileName);

%% show the origenal image and its gray value
 figure, imshow(Im), title('The RGB Image');
 figure, imshow(rgb2gray(Im)), title('The gray value Image');

%% show different channel of RGB

 figure, imshow(Im(:,:,1)), title('The Red channel');
 figure, imshow(Im(:,:,2)), title('The Green channel');
 figure, imshow(Im(:,:,3)), title('The Blue channel');

%% change the colorspace to Lab 

[l a b] = RGB2Lab(Im);

 figure, imagesc(l), colormap(gray), title('The l channel in CIELAB colorspace');
 figure, imshow(a), title('The a channel in CIELAB colorspace');
 figure, imshow(b), title('The b channel in CIELAB colorspace');

%% change the colorspace to HSV

hsv_image = rgb2hsv(Im);

 figure, imshow(hsv_image(:,:,1)), title('The H plane in HSV colorspace');
 figure, imshow(hsv_image(:,:,2)), title('The S plane in HSV colorspace');
 figure, imshow(hsv_image(:,:,3)), title('The V plane in HSV colorspace');

%% apply fast radial symmetry transfrom to extract cells
hsv_s255 = uint8( hsv_image(:,:,2) * 255);

%smooth the image 
gaussK = fspecial('gaussian',[7 7],2);
smoothhS255v = imfilter(hsv_s255, gaussK,'same');

[S, So] = fastradial(smoothhS255v, [3 4 5], 1, .2);
figure 
subplot(2,1,1), imshow(S), title('a) FRST of S pplane in HSV');
% since all cells are bright in S, then the positive value in FRST indication 
% nuclei/cell
bin_s = im2bw(S,graythresh(S));
subplot(2,1,2), imshow(bin_s), title('binary image of a');
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
Blue= Im(:,:,3);
Blue(BinarySOpenByRec) = 255;
figure, imshow(Blue)

% figure
% imshow(Im), hold('on');
% contour(BinarySOpenByRec)

%% smoothing with gaussian kernel
gaussK = fspecial('gaussian',[7 7],2);
Ig = imfilter(hsv_s255, gaussK,'same');
figure, imshow(Ig), title('smoothed saturation with Gaussian');
%%
SE2 = strel('disk', 2); % structure element for erosion to dis connect
ErodedAChannel = imerode(a, SE2);
figure, imshow(ErodedAChannel), 
title('The eroded a channel with disk SE of radius 3');

doubleEroded = imerode(ErodedAChannel, SE2);
RecbyOpen = imreconstruct(doubleEroded, ErodedAChannel);
figure, imshow(RecbyOpen), 
title('Reconstruction by opening');
BW = im2bw(RecbyOpen, graythresh(RecbyOpen));
figure, imshow(BW), title('The binary image');


%% new colorspcae 

maxI = max(Im, [], 3);
meanI = mean(Im,3);

%% try to get only germ mass by using color plane that the ger mmass
% are shown more discriminatively (like red in RGB, l in lab and v hsv).
% note: that grem mass are dark region in these channals

% l plane in lab seems to be the baest one 
% folder = '/home/rahim/work/Dropbox/images/GATA-4/';
% files = dir([ folder '*.TIF']);
% for i = 1: numel(files)
%     imgFile = [folder files(i).name];
%     Im = imread(imgFile);
    [l a b] = RGB2Lab(cartoonI2); % Im 
    Red = percentileStretch(l, [0, 1]);% Im(:,:,1);
    figure(100412), imshow(Red), title('The chosen channel');
    %Red = uint8(gaussfilt(Red, 1));
    
    % anisRed = percentileStretch(anisodiff(Red, 100, 40, .25, 2), [0,1]);
    % inverse of binary of red channel
    BWRed = not(im2bw(Red, graythresh(Red)));  %not
    figure(100112), imshow(BWRed);
    title('gray threshold the inverse of red channal');
    
    lblRed =label(BWRed,1, 100,0);
    
    % eroded7 = erosion(lblRed>0 , 7); best config for DIPimage
    
    lblRed = im2mat(lblRed >1);
    SE4= strel('disk', 4);
    erodedRed =imerode(lblRed, SE4);
    RecRed = imreconstruct(erodedRed, lblRed);
    figure(100212), imshow(RecRed, []);
    title('the remain region of red channal after reconstruction');
    
    % calculate distance transfrom to connect regions
    [dist nearstlbl] = bwdist(RecRed);
    % only connect region that their distance lower than 10 pixels
%     figure(3), imshow(dist < 10), title('threshold of distance transfrom at 10');
    figure(100312), imshow(Red);
    RGBDist = uint8(cat(3, zeros(size(Red)),(dist<10)*255 ,zeros(size(Red))));
    hold on
    distIm = imshow(RGBDist); title('regions on top of the image');
    set(distIm,'AlphaData',0.3);

%     pause
% end



