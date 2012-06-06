function out = val2obj(Im, val)
% function out = val2obj(Im, val)
% Map value to the labeled region of an image.
% Input,
%   Im, a labeled image.
%   val, an 1D array that each element of the array has a value for its 
%   coresponding region.
% Output,
%   out, an image with that same size as Im, but hte regoins are labeled
%   with their corresponding value.
% 
% Author: A.Rahim Kadkhodamohammadi (r.k.mohammadi@gmail.com)
% 01 June 2012 CBA, Uppsala University
%--------------------------------------------------------------

sz = size(Im);
out= zeros(sz(1:2));
for row = 1 : sz(1)
    for col = 1 : sz(2)
        if Im(row,col)
            out(row,col) = val(Im(row,col));
        end
    end
end

end