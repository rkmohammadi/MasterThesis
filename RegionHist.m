function regions = RegionHist(I, stats, quanNum)
% function regions = RegionHist(I, stats, quanNum)
% claculate the histogrma of each region that its pixel coordinate is 
% specified in stats.PixelList. the histograms are calculated per each 
% channel of input image (I).
% inputs: 
%       I, input image or 3D matrix of channels that we want calculate 
%          histogram there. 
%       stats, structure array that conntain the pixel list of eahc region
%              it could be cmoputed by regionprops
%       quanNum, the number of bin in histogram i.e quantization level,
%                default value is 255.
% output:
%       regions, array of structure that contain the histogram of each
%                region. E.x the histogram of plane i of input I in region
%                r can be accessed via regions(r).hist(i,:)
% 
% author: A.Rahim Kadkhodamohammadi r.k.mohammadi@gmail.com
% February 14, 2012
%--------------------------------------------------------------

if nargin < 2 || nargin > 3
    error('Please check function description for the number of input');
elseif nargin == 2
    quanNum = 255;
end

sz =size(I);
bins = zeros(quanNum, sz(3));
for i=1 : sz(3)
    bins(:,i) = linspace(min(min(I(:,:,i))),max(max(I(:,:,i))), quanNum);
end
numRegion = length(stats);
s = struct('hist', zeros(sz(3), quanNum));
regions = repmat(s, numRegion, 1);
for dim =1 : sz(3) % the number of plane in input image
    for i=1 : numRegion
        PL = stats(i).PixelList;
        rows = size(PL,1); % the number of pixels in the region
        pixels = zeros(rows,1);
        % extract region's pixel
        for j=1 : rows
            pixels(j) = I(PL(j,2), PL(j,1), dim);
        end
        %sote the normilized histogram of each region per channel
        regions(i).hist(dim,:) = hist(pixels, bins(:,dim))./rows; 
    end
end