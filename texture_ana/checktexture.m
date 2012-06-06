%% doing some experiment with LBP
% fisrt select some region and then apply lbp on those regions

fileName = '../../images/set_1/M1x10ex1.tif';
Im = imread(fileName);
Im = imresize(Im,0.2,'lanczos2');
% selected channel 
I = Im(:,:,3);

% show figure to select for different region
choice = input('load default mask (y/n)?', 's');
if choice == 'y'
    load('m02.mat');
else
    figure(1011);
    %msgbox('Please select Lumen');
    masks(1).BW = roipoly(I);
    
    figure(1011);
    %msgbox('Please select germ cell region');
    masks(2).BW = roipoly(I);
    
    figure(1011);
    %msgbox('Please select sertoli');
    masks(3).BW = roipoly(I);
    
    figure(1011);
    %msgbox('Please select interstitial');
    masks(4).BW = roipoly(I);
    close(1011)
end

% extract bounding box
numMask = numel(masks);
for i =1 : numMask
    stats = regionprops(masks(i).BW, 'BoundingBox');
    masks(i).boundingBox = round(stats.BoundingBox);
end

%mapping for LBP
R = 2; % 1
N = 8; % 8
mapping=getmapping(N,'ri');

% texture in R-G-B
double_Im = double(Im);
figure('name', 'RGB space');
for i = 1 : 3
     for j=1 : numMask
         currentROI = double_Im(:,:,i) .* masks(j).BW;
         currentROI = currentROI(masks(j).boundingBox(2): ...
             (masks(j).boundingBox(2)+masks(j).boundingBox(4)), ...
             masks(j).boundingBox(1): ...
             (masks(j).boundingBox(1)+masks(j).boundingBox(3)));
%          nH=lbp(currentROI,R,N,mapping,'nh');
         nH=lbp_masked(currentROI,currentROI>0,R,N,mapping,'nh');
           
%          figure, subplot(2,1,1), imshow(double_Im(:,:,i) .* masks(j).BW);
%          tmp = double_Im(:, : ,i);
%          tmp(masks(j).boundingBox(2): ...
%              (masks(j).boundingBox(2)+masks(j).boundingBox(4)), ...
%              masks(j).boundingBox(1): ...
%              (masks(j).boundingBox(1)+masks(j).boundingBox(3)),i) =255;
%          subplot(2,1,2), imshow(tmp);
         
         subplot(3,numMask, (i-1)*numMask + j), stem(nH);
     end
end

% H-S-V
hsv_image = rgb2hsv(Im);
figure('name', 'HSV space');
for i = 1 : 3
     for j=1 : numMask
         currentROI = hsv_image(:,:,i) .* masks(j).BW;
         currentROI = currentROI(masks(j).boundingBox(2): ...
             (masks(j).boundingBox(2)+masks(j).boundingBox(4)), ...
             masks(j).boundingBox(1): ...
             (masks(j).boundingBox(1)+masks(j).boundingBox(3)));
         
         nH=lbp_masked(currentROI,currentROI>0,R,N,mapping,'nh');
         subplot(3,numMask, (i-1)*numMask + j), stem(nH);
     end
end

% L-a-b
lab_image = RGB2Lab(Im);

figure('name', 'Lab space');
for i = 1 : 3
     for j=1 : numMask
         currentROI = lab_image(:,:,i) .* masks(j).BW;
         currentROI = currentROI(masks(j).boundingBox(2): ...
             (masks(j).boundingBox(2)+masks(j).boundingBox(4)), ...
             masks(j).boundingBox(1): ...
             (masks(j).boundingBox(1)+masks(j).boundingBox(3)));
         
         nH=lbp_masked(currentROI,currentROI>0,R,N,mapping,'nh');
         subplot(3,numMask, (i-1)*numMask + j), stem(nH);
     end
end

% I = ROI .* double(I);
% stats = regionprops(ROI, 'BoundingBox');
% boundingBox = stats.BoundingBox;
% I = I(boundingBox(2): (boundingBox(2)+boundingBox(4)), ...
%     boundingBox(1): (boundingBox(1)+boundingBox(3)));
% 
% mapping=getmapping(8,'ri');
% H1=lbp(I,1,8,mapping,'nh'); %LBP histogram in (8,1) neighborhood
% %using uniform patterns
% figure
% subplot(2,1,1),stem(H1);
% 
% H2=lbp(I);
% subplot(2,1,2),stem(H2);
