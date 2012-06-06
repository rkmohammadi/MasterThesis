function [coords, grayVals] = get_coords(points, im, method)
% function [coords, grayVals] = get_coords(points, im)
%  this function get coordinate of two points and return the coordinate of
%  pixels that connect this two points. If you call this method with second
%  argument it will return an array of gray value corresponding to those
%  pixels. By third argument you can determine the method of interpolating
%  value of pixels, defalut value is 'nearest'.
% input:
%       points, coordinate of two points [x1 y1;x2 y2]
%       im(optional), image to find the gray values of pixels in between points
%       method(optional), 'linear', 'spline', 'cubic', 'nearest'
% output:
%       coords, coordinate of pixels in between points
%       grayVals, gray values of coords
% 
% A.Rahim Kadkhodamohammadi (r.k.mohammadi@gmail.com)
% February 7 /2012
%--------------------------------------------------------------------------

if nargin <3
    method = 'nearest';
end

N = size(points,1);
if N ~= 2
   error('You need to pass exactly two points.')
end

len = norm(points(2,:)-points(1,:));
k = floor(len);
coords = repmat(points(1,:),[k,1])+(0:k-1)'*((points(2,:)-points(1,:))/len);
 
% return the gray value or color valur of the pixels between points
if nargin > 1 
    ch = size(im,3);
    if ch ==1
        grayVals = get_subpixel(squeeze(im),coords, method);
    else
        grayVals = zeros(k,ch);
        for ii =1 : ch
            grayVals(:,ii) = get_subpixel(squeeze(im(:,:,ii)),coords, method);
        end
    end
end

end