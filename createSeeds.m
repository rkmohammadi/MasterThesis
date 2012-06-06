function createSeeds(img)
% function createSeeds(img)
%  A function  that get a file name and show the file to chose regions and
%  the selected regions are stored in an image file with the same name as
%  the input file in PNG format.
% 
% Author: A.Rahim Kadkhodamohammadi (r.k.mohammadi@gmail.com)
% 01 June 2012 CBA, Uppsala University
%--------------------------------------------------------------

if nargin == 0
    img = '/home/rahim/work/Dropbox/images/GATA-4/M1G40024.TIF';
elseif nargin > 1
    error('this function only get one string as image filename');
end

I = imread(img);
Rseeds = I(:,:,1);
Gseeds = I(:,:,2);
Bseeds = I(:,:,3);

inc = 'y';
while inc == 'y'
    BW = roipoly(I);
    color= input ('input color for the selected region: ');
    Rseeds(BW) = color(1);
    Gseeds(BW) = color(2);
    Bseeds(BW) = color(3);
    inc = input('do you want select more region (y to contiue)? ', 's');
end
seeds = cat(3, Rseeds, Gseeds, Bseeds);
imwrite(seeds, [img(1:end-3) 'png'], 'png');
end