function out = combineGMLumen(GM, PLR, numL)
% function out = combineGMLumen(GM, PLR, numL)
%  this function get germ-mass and potential lumen regions(PLR) and look at
%  their intersection. We expand the PLR by distThre (defalut is 50) in 
%  each direction. We also ignore germ-mass smaller than sizeThre (5000).
%  The other parameter are defined on my report section "Detect Lumen".
%  Input:
%        GM, a binary image contains the germ-mass regions;
%        PLR, a binary image contains Potential Lumen Region and its size
%           is equal to GM;
%        numL, the number of separate region in PLR, i.e. the number of
%        potential lumen regions.
% 
%  Output:
%       out, an binary image with the same size as GM which contains the
%       accepted lumen.
% 
% Author: A.Rahim Kadkhodamohammadi (r.k.mohammadi@gmail.com)
% 01 June 2012 CBA, Uppsala University
%--------------------------------------------------------------

if nargin ~= 3
    error('Please check function help for appropriate inputs!');
end

% the paremeter of the algorithm that can be tuned here.
distThre = 50;
sizeThre = 5000;

% this comstant are based on section detect Lumen on my report eq. 4.6 and
% 4.7
beta1 = 0.85; % when we have only one germ-mass around the PLR
beta = 0.6; 
gamma =  0.8;


sz = size(GM);
out = zeros(sz(1:2));

[GMlbl numGM]= bwlabel(GM);
stats = regionprops(GMlbl, 'Perimeter', 'Area');
% distance threshold
for l = 1 : numL
    currentL = (PLR == l);
    distL = bwdist(currentL);
    growDist= (distL < distThre);
    boundary = bwboundaries(growDist, 'noholes');
    lenb= size(boundary{1},1);
    % remove pixeles that place in the boarder of image
    boarderP = length(find(boundary{1}(:,1)==1));
    boarderP =  boarderP + length(find(boundary{1}(:,1)==sz(1)));
    boarderP =  boarderP + length(find(boundary{1}(:,2)==1));
    boarderP =  boarderP + length(find(boundary{1}(:,2)==sz(2)));
    lenb = lenb - boarderP;
    
    touchedR = and(growDist,GMlbl);
    
    [boundaryTouchedR lblTR]= bwboundaries(touchedR, 'noholes');
    if numel(boundaryTouchedR) == 1
        if size(boundaryTouchedR{1},1)/2 > beta1*lenb % .9
            tt = and(growDist, ~touchedR);
            tt = bwareaopen(tt, 5000);
            out = or(out, tt);
            clear tt
        end
    else % there more than one neighbor germ mass
        lenbTR = 0;
        for ii = 1 : length(boundaryTouchedR)
            if size(boundaryTouchedR{ii},1) > 50
                lenbTR = lenbTR + size(boundaryTouchedR{ii},1)/2-10;
            end
        end
        if lenbTR > beta*lenb
            combine = 1;
            GM2Boundary =[];
            for ii = 1 : length(boundaryTouchedR)
                [~,~, GMlblId] = find(GMlbl(lblTR == ii), 1, 'first');
                ind = [];
                if numel(GM2Boundary)
                    ind = find(cat(1,GM2Boundary.lbl) == GMlblId);
                end
                if isempty(ind)
                    GM2Boundary = [GM2Boundary, struct('lbl',GMlblId,'bind',ii)];
                else
                    GM2Boundary(ind).bind = [GM2Boundary(ind).bind , ii];
                end
            end
            for ii = 1 : numel(GM2Boundary)
                GId = GM2Boundary(ii).lbl;
                szBoundary = 0;
                bInd = GM2Boundary(ii).bind;
                for jj = 1 : numel(bInd)
                    szBoundary = szBoundary + size(boundaryTouchedR{bInd(jj)},1)/2;
                end
                if stats(GId).Area > sizeThre && ...
                        szBoundary >= (gamma * stats(GId).Perimeter)
                    combine =0;
                    break;
                end
            end
            if combine
                tt = and(growDist, ~touchedR);
                tt = bwareaopen(tt, 5000);
                out = or(out, tt);
                clear tt
            end
        end 
    end
    
end
end