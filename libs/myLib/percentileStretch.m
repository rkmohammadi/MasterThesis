function IStretched = percentileStretch(I, percentiles)
% percentile stretch default value for percentile is [.05 .95]
% input image could be gray or color image and the type output 
% image is uint8. 
%
% Author: A.Rahim Kadkhodamohammadi r.k.mohammadi@gmail.com
% 25/01/2012 CBA, Uppsala University
%--------------------------------------------------------------

if nargin == 1 
    percentiles = [.05, .95];
end

if size(I,3) == 1
    % number of pixel
    N = length(I(:));
    IStretched = PerStretch(I, percentiles, N);
else
    N = numel(I(:,:,1));
    for i =1 : size(I,3)
        IStretched(:,:,i) = PerStretch(I(:,:,i),percentiles, N);
    end
end

end
function IStretched = PerStretch(I, percentiles, N)
I = double(I);

% bring the image to [0 255] interval 
low = min(min(I));
high = max(max(I));
I = round((I -low) .* (255/(high - low)));
H = histc(I(:), 0:255);

cutoffs = round(percentiles .* N);

cumHist = cumsum(H);

a = find(cumHist >= cutoffs(1));
a = a(1) -1;

b = find(cumHist >= cutoffs(2));
b = b(1) -1;


R = b-a;
scale = 255/R;
IStretched = round((I-a) .* scale);
IStretched = uint8(IStretched);
end
