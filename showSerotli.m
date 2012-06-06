function showSerotli(cond, CloseGM, imgOut, h, sertoliLbl, GMandSertoli, en_im, posAndNeg) 
% function showSerotli(cond, CloseGM, imgOut, h, sertoliLbl, GMandSertoli, en_im, posAndNeg) 
% A function to show the Sertoli cells and germ-mass based on the
% condition. More information is given in section 5.2 and 5.3.
% Input: 
%   cond, bolean if one check angle just if zero consider ang for connection
%   that crosses one or two edges in boundary map.
%   CloseGM, the neighbor infomation for the Sertoli cells computed by
%   findCloseGM.
%   imgOut, a filename to store the output result
%   h, handle to a figure to show the result over the figure, if you want
%   to show the image in new figure you need to provide the rest of the
%   arguments
%   sertolilbl, a label image of sertoli cells
%   GMandSertoli, a label image of germ-mass and Sertoli cells.
%   en_im, a color image 
%   posAndNeg, the structure that have the neighor of each Sertoli cells
%   per each side of connecting plane
% 
%  See also findClosestGM, combineGMLumen, findNeighbors.
% 
% Author: A.Rahim Kadkhodamohammadi (r.k.mohammadi@gmail.com)
% 01 June 2012 CBA, Uppsala University
%--------------------------------------------------------------

if nargin > 5
    [~, h] = overlayCluster(en_im, GMandSertoli, [], 1); 
end
hold on
% sertoliLbl = bwlabel(sertoli,4);
stats = regionprops(sertoliLbl, 'Centroid');
num = numel(stats);

for i = 1 : num
    NS = CloseGM{i};
    for j =1 : size((NS{2}),1)
        if isempty(NS{3}{j}) %&& NS{2}(j,5)< 3 %&& NS{2}(j,5) == 0 
            
           if NS{2}(j,5) < 2 || cond
            line([stats(i).Centroid(1), NS{2}(j,4)], ...
                [stats(i).Centroid(2), NS{2}(j,3)]);
            else
                line([stats(i).Centroid(1), NS{2}(j,4)], ...
                [stats(i).Centroid(2), NS{2}(j,3)], 'color', 'r');
            end
            if ~isempty(NS{1}) && checkCondition(i,stats,posAndNeg,j, NS{1}, NS{2}(j,1),NS{2}(j,2) ,NS{2}(j,5), cond)
                plot(stats(i).Centroid(1), stats(i).Centroid(2), ...
                    'Marker','p','Color',[.88 0.75 .99],'MarkerSize',20);
            end
            break
        end
    end
end
saveas(h, imgOut);

end



function out = checkCondition(currentS, stats,posAndNeg, ind, NSId, GMId, GMDist, numBorders ,cond)
out = 0;
try
    neighborInfo = posAndNeg{currentS}(ind);
catch err
    disp(err);
    disp(['cuurevtS:' int2str(currentS)]);
    return
end
if isempty(neighborInfo.neg) || isempty(neighborInfo.pos)
    return
end
negId = NSId(neighborInfo.neg);
posId = NSId(neighborInfo.pos);
negcntr = cat(1,stats(negId).Centroid);
poscntr = cat(1,stats(posId).Centroid);

distNegs = sqrt(sum((negcntr-repmat(stats(currentS).Centroid,numel(negId),1)).^2,2));
distPoss = sqrt(sum((poscntr-repmat(stats(currentS).Centroid,numel(posId),1)).^2,2));

closeNeg = closest(distNegs,GMId, neighborInfo.negGMlbl);
closePos = closest(distPoss, GMId, neighborInfo.posGMlbl);
if closeNeg == 0 || closePos == 0
    return
end
ang = calculateAng(negcntr(closeNeg,:), stats(currentS).Centroid, poscntr(closePos,:));
distNeg = neighborInfo.negGMDist(closeNeg);
distPos = neighborInfo.posGMDist(closePos);
if abs(ang)<1 && (cond ||numBorders==1) %||(numBorders==1 && (mean([distNeg , distPos]) < 1.2 *GMDist) )%  ang>0 && ang<.8
    out = 1;
end

end

function ind = closest(dist,GMId, GMlbl)
[sortDist, sortId] = sort(dist);
ind  = 0;
for i =1 : length(sortDist)
    if GMlbl(sortId(i)) == GMId
        ind = sortId(i);
        return
    end
end
end

function ang = calculateAng(c1, c2, c3)
a= norm(c1-c2);
b= norm(c2 -c3);
c= norm(c1-c3);
ang = (a^2 + b^2 - c^2)/(a*b);
end