function  [clusterdRgn h] = overlayCluster(I,lbl_im,lbl2cluster, numCluster,...
    imgFile, germCluster)
% function  clusterdRgn = overlayCluster(I,lbl_im,lbl2cluster, numCluster,...
%     imgFile, germCluster)
%  show cluster on top of the input image. if the number of cluster
%  (numCluster) equals one, the lbl_im will be shown on top of I. 
% input:
%       I, image that want ot show cluters on top of it. it could be color
%       or gray value image. Now this function works well for color image.
%       lbl_im, labeled image s.t. if sz = size(I), size(lb_im) == sz(1:2).
%       lbl2cluster, 1D array of cluster indicies that lbl2cluster(i) is
%       the cluster number of region with label i.
%       numCluster, scaler 1*1 that is the number of cluster (optional). If
%       this value set to  one, lbl2cluster will be ignored and lbl_im will
%       beshowed on top of I.
%       imgFile, string that contain the image file name (optional) or cell
%       array with two element ('save', filename) to save dispalaied image
%       with the filename.
%       germCluster, cluster number of germ mass.
% ouputp:
%       clusterdRgn, matrix with the same size as lbl_im that cluster number
%       has assigned regions.
%       h, figure handle for the last image
% 
% A.Rahim Kadkhodamohammadi (r.k.mohammadi@gmail.com)
% March 13 /2012, CBA Uppsala University
%--------------------------------------------------------------------------
storeFlag = 0;

if nargin < 3
    error('This function needs at least three inputs !!!');
elseif nargin < 4
    numCluster = max(lbl2cluster(:));
    imgFile = '';
elseif nargin < 5
    imgFile = '';
elseif nargin >= 5
    if ischar(imgFile)
        imgFile = ['(' imgFile ')'];
    else
        if iscellstr(imgFile) && numel(imgFile) == 2 && strcmpi(imgFile{1}, 'save')
            storeFlag = 1;
            imgFile = imgFile{2};
        end
    end
end

if nargin < 6
    germCluster = -1;
end

sz = size(I);
if numCluster == 1
    clusterdRgn = (lbl_im > 0);
else
    clusterdRgn= zeros(sz(1:2));
    for row = 1 : sz(1)
        for col = 1 : sz(2)
            if lbl_im(row,col)
                clusterdRgn(row,col) = lbl2cluster(lbl_im(row,col));
            end
        end
    end
end

for c = 1: numCluster
    h = figure; 
    imshow(I);
    RGBDist = uint8(cat(3, zeros(sz(1:2)),(clusterdRgn == c)*255, ...
        zeros(sz(1:2))));
    
    hold on
    distIm = imshow(RGBDist);
    if c == germCluster
        title([imgFile ' Cluster ' 'Germ Mass']);
    else
        title([imgFile ' Cluster ' int2str(c) ' out of ' int2str(numCluster)]);
    end
    set(distIm,'AlphaData',0.3);
    hold off
    if storeFlag
        saveas(h, [imgFile(1:end-4) 'c' num2str(c) imgFile(end-3:end)]);
    end
end
end