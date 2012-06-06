function  [totalGfe, totalNonGfe]= showfeaturesofSeeds(verbose)
% function  [totalGfe, totalNonGfe]= showfeaturesofSeeds(verbose)
% This funcion read following image and their stored labled cells to
% compute some feature and show them. If call this function and set verbose
% to one it will the labled image over the original one.
%
% Author: A.Rahim Kadkhodamohammadi (r.k.mohammadi@gmail.com)
% 01 June 2012 CBA, Uppsala University
%--------------------------------------------------------------

if nargin <1 
    verbose = 0;
end
% manually define class of some region
% the file name are stored invariable images

images = {'../images/GATA-4/M3G4001.TIF'; ...
    '../images/GATA-4/M3G4004.TIF'; ...
    '../images/GATA-4/M3G4006.TIF'; ...
  %  '../images/GATA-4/M1G40019.TIF'; ...
    '../images/GATA-4/M1G40022.TIF'; ...
    '../images/GATA-4/M1G40036.TIF'};

totalGfe = [];
totalNonGfe = [];
for i = 1 : numel(images)
    imgFile = images{i};
    load([imgFile(1:end-3) 'mat']);
    germS = class2Region(1:end-1);
    nonGermS = [class1Region(1:end-1); class4Region(1:end-1)];

    if  verbose
       showseeds(germS, lbl_bin_FRST, en_im, imgFile);
       showseeds(nonGermS, lbl_bin_FRST, en_im, imgFile);
    end 

    measurmentIDs = {'Size', 'Perimeter', 'Mean', 'StdDev', ...
        'Skewness', 'MaxVal', 'MinVal', 'Center', 'Inertia','Center'};
    lab = RGB2Lab(en_im);
    stretchedL = percentileStretch(lab(:,:,1), [0,1]);
    diplbl = label(lbl_bin_FRST>0,1,20,0);
    msr = measure(diplbl, mat2im(stretchedL),measurmentIDs, [],1);
    
    Gfe = [msr(germS).Skewness; msr(germS).Mean; msr(germS).StdDev; msr(germS).MinVal];
    nGfe = [msr(nonGermS).Skewness; msr(nonGermS).Mean; msr(nonGermS).StdDev; msr(nonGermS).MinVal];
    
    totalGfe = [ totalGfe Gfe];
    totalNonGfe = [ totalNonGfe nGfe];
end
figure, plot3(totalGfe(4,:), totalGfe(2,:), totalGfe(3,:),'.');
hold on, plot3(totalNonGfe(4,:), totalNonGfe(2,:), totalNonGfe(3,:),'x');
end


function showseeds(SeedRegion, lbl_bin_FRST, en_im, imgFile)
figure,
imshow(en_im);
hold on;
mask = zeros(size(lbl_bin_FRST));
for j =1 : length(SeedRegion)
    if SeedRegion(j) ~= 0
        mask = mask + (lbl_bin_FRST== SeedRegion(j));
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