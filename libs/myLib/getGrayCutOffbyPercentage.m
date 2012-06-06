function level = getGrayCutOffbyPercentage(img, percentage)
% function level = getGrayCutOffbyPercentage(img, percentage)
%  this will get and a gray value image and return the a gray value that greater
%  than or equal to percentage of gray value of the image.
% input:
%       img, gray value image
%       percentage, a scaler value between [0 1) to determine the cut-off
%       value, default is 0.9.
% ouput:
%       level, gery level that greater than or equal percentage of image.
% 
% Author: A.Rahim Kadkhodamohammadi r.k.mohammadi@gmail.com
% 23/03/2012 CBA, Uppsala University
%--------------------------------------------------------------

if nargin == 1
    percentage = .9;
end
I = round(img(:)*100);
N = length(I);

cutoff = percentage * N;
bins = 0:max(I);
H = histc(I, bins);
cumHist = cumsum(H);
a = find(cumHist >= cutoff, 1,'first');
level = bins(a-1)/100;
end
