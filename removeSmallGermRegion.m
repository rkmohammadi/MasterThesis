function out = removeSmallGermRegion(grayImg, germMass, cutOffs)
% function out = removeSmallGermRegion(grayImg, germMass, cutOffs)
% A function to tune the germ-mass regions according to their size and
% intensity.
%  input:
%   grayIm, an gray level image.
%   germMass, a binary image of germ-mass regions.
%   cutOffs, an array of threshold value such as
%   1*4 array; first value maximum size of small region to remove, 
%   second value the maximum size fo regions that should be checked 
%   third value percentage of high gray value that should be in the region
%   to remove, fourth precentage to determine the high value gray level;
%   for example: [800, 7000, 0.6, 0.9].
%  Note: percentage should be in [0,1)
%  Output:
%   out, a binary image of the adjusted germ-mass regions.
% 
% Author: A.Rahim Kadkhodamohammadi (r.k.mohammadi@gmail.com)
% 01 June 2012 CBA, Uppsala University
%--------------------------------------------------------------

sz = size(germMass);
level = getGrayCutOffbyPercentage(grayImg, cutOffs(4));
minSize = cutOffs(1);% 800
% germ mass region greater that threshold that are acceptable
GM1 = bwareaopen(germMass, minSize);  
acceptableGM = bwareaopen(germMass, cutOffs(2));
roi = and(GM1, ~acceptableGM);
[lblRoi num] = bwlabel(roi);
stats = regionprops(lblRoi, 'Area');
s1 = cat(2,stats.Area)*cutOffs(3);

lowGrayValInd = grayImg < level;
tmpLblRoi = lblRoi;
tmpLblRoi(lowGrayValInd) = 0;
% set the i pixel for each label to haveat least one pixel for each of them
% to correctly compute area.
tmpLblRoi(1:num,1) = 1:num; 

tmpStats = regionprops(tmpLblRoi, 'Area');
s2 = cat(2, tmpStats.Area);
lblPassedRegion = find(s2 < s1);
numPassedR = numel(lblPassedRegion);
if numPassedR < num
passedR = zeros(sz(1:2));
for i = 1 : numel(lblPassedRegion)
    passedR = or(passedR, lblRoi==lblPassedRegion(i));
end

out = or(passedR,acceptableGM);
else
    out = GM1;
end