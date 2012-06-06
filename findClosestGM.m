function [closeGM sertoliLbl centers GMandSertoli posNeg] = findClosestGM(sertoli, connectedGerm, ND, Im, thre)
% function [closeGM sertoliLbl centers GMandSertoli posNeg] ...
%               = findClosestGM(sertoli, connectedGerm, ND, Im, thre)
% This will get the Sertoli cells and germ-mass to build a data structure
% that hold infomation such as closest germ-mass, neighboring sertoli cells
% etc. More information can be found in 5.1 "Assign non-germ cells to germ-mass"
%  Input;
%   sertoli, a binary image contains the Sertoli cells.
%   connectedGerm, a binary image of germ-mass.
%   ND, an integer, maximum neighborhood distance.
%   Im (optional), an image with the same size as sertoli.
%   thre(optional needed if Im is given), a threshold value for
%   thresholding the Im to count the number crosses when connecting a cell
%   to an germ-mass.
%  Ouput:
%   closeGM, is a cell array that has entry for each Sertoli cell and each
%   entry is array of 3 cell; First, a column array stores the label indices
%   of neighboring Sertoli cells in the specified neighborhood distance;
%   Second 2D array that each row has germ-mass label, distance to closest 
%   point in the border of germ-mass, row and column indices of that point,
%   and the number of white pixels thresholded Im if it has given; I need 
%   to mention that this matrix is sorted based on the cell's distance to 
%   the germ-mass;Third one has tha same number of row as the second one 
%   and store the indices of srtoli celles that are in-between the current
%   sertoli and germ mass;
%   Note: it is good to replace the cell array with structure array which
%   is more convinient.
%   sertoliLbl, an image contains the labeled Sertoli cells.
%   centers, a 2D array that each row has the centroid of the corresponding
%   Sertoli cells.
%   GMandSertoli, a label matrix that the Sertoli cells are labeled as
%   sertoliLbl and the Germ-mass regions are labled with negetive numbers.
%   posNeg, an structure that have two field 'pos' and 'neg' that contains
%   information for all neighbors of the Sertoli cells; the connection
%   between the sertoli cell and the Germ-mass splite the plane into two
%   negetive and positive plane that the neighbors of the cell in each
%   plane are stored in corresponding feild.
% 
% Author: A.Rahim Kadkhodamohammadi (r.k.mohammadi@gmail.com)
% 01 June 2012 CBA, Uppsala University
%--------------------------------------------------------------

sertoliLbl = bwlabel(sertoli,4);
germBoundary = imdilate(connectedGerm,strel('disk',1)) - connectedGerm;
germBoundaryLbl = bwlabel(germBoundary,8);
stats = regionprops(sertoliLbl, 'Centroid');
num = numel(stats);
% cell of 2D that the first column store a list of direct germ mass
% neighbors and their distance. The second column store the list of germ
% Mass that connected to the current sertoli cell via another sertoli
% cells. So, it includes germ mass label, distance and connecting sertoli.
closeGM = cell(num,1);
posNeg = cell(num,1);

GMandSertoli = sertoliLbl - germBoundaryLbl;
[row, col] = size(GMandSertoli);
centers = round(cat(1,stats.Centroid));
switch nargin
    case 3
        for i = 1 : num
            cntr = centers(i,:);
            NS = neighboringSertoli(i, cntr, row,col,GMandSertoli,ND);
            if numel(NS{1}) > 0
                directNS = directNeighbor(i, centers(i,:), NS{1}, centers, ...
                    sertoliLbl, row, col);
                
                NS{1} = directNS;
            end
            closeGM{i} = NS;
        end
        for i = 1 : num
            cntr = centers(i,:);
            NS = closeGM{i};
            if numel(NS{1}) > 0
                [nGermMass ~] = size(NS{2});
                if nGermMass > 0
                    posNeg{i} = separateNeighbor(nGermMass, NS, cntr, centers, closeGM);
                end
            end
        end
    case 5
        % give each sertoli germ mass a weight which is the portion of pixels
        % above threshold per length of connecting line
        BW = Im > thre;
        for i = 1 : num
            cntr = centers(i,:);
            NS = neighboringSertoliWeight(i, cntr, row,col,GMandSertoli,ND, BW);
            if numel(NS{1}) > 0
                directNS = directNeighbor(i, centers(i,:), NS{1}, centers, ...
                    sertoliLbl, row, col);
                
                NS{1} = directNS;
            end
            closeGM{i} = NS;
        end
        for i = 1 : num
            cntr = centers(i,:);
            NS = closeGM{i};
            if numel(NS{1}) > 0
                [nGermMass ~] = size(NS{2});
                if nGermMass > 0
                    posNeg{i} = separateNeighbor(nGermMass, NS, cntr, centers, closeGM);
                end
            end
        end
        
    otherwise
        error('Please check the function for correct number of input');
end

end

function NS = neighboringSertoli(currentR, centroid, row,col,lbl,ND)
currentBox = lbl(max(1,centroid(2)-ND): min(centroid(2)+ND, row), ...
    max(1,centroid(1)-ND): min(centroid(1)+ND, col));
[~, ~, v] = find(currentBox);
NS = {[], [], []};
v = sort(v);
for j = 1 : numel(v)-1
    if (v(j) ~= currentR && v(j) ~= v(j+1))
        if v(j) > 0
            NS{1} = [NS{1}, v(j)];
        else
            NS{2} = [NS{2}, v(j)];
        end
    end
end

if (numel(v) > 0 && v(end) ~= currentR )
    if v(end) > 0
        NS{1} = [NS{1}, v(end)];
    else
        NS{2} = [NS{2}, v(end)];
    end
end
numGM = numel(NS{2});
if numGM > 0
    tmpNS = zeros(numGM,4);
    pixelPos = zeros(numGM,2);
    cntrInROI = min(centroid, [1 1]+ND);
    for l = 1 : numGM
        [r,c] = find(currentBox == NS{2}(l));
        D = sqrt((r-cntrInROI(2)).^2 + (c-cntrInROI(1)).^2);
        [val, ind] =min(D);
        tmpNS(l,1) = - NS{2}(l); % germ mass label
        tmpNS(l,2) = val; % ditance of centriod and nearest pixel in boundary
        % first coordinate of nearest pixel or row
        tmpNS(l,3) = r(ind)+ max(1,centroid(2)-ND)- 1;
        % second coordinate of nearest pixel or column
        tmpNS(l,4) = c(ind)+max(1,centroid(1)-ND) - 1;
        pixelPos(l, :) = [c(ind) r(ind)];
    end
    % sort neighboring germ mass based on distance
    [~, sortTmpNS] = sort(tmpNS(:,2)); 
    NS{2} = tmpNS(sortTmpNS, : );
    pixelPos = pixelPos(sortTmpNS, :);
    NS{3} = checkBridgeRegion(currentR, cntrInROI, NS{2},pixelPos, currentBox);
end
end


function NS = neighboringSertoliWeight(currentR, centroid, row,col,lbl,ND, BW)
currentBox = lbl(max(1,centroid(2)-ND): min(centroid(2)+ND, row), ...
    max(1,centroid(1)-ND): min(centroid(1)+ND, col));
[~, ~, v] = find(currentBox);
NS = {[], [], []};
v = sort(v);
for j = 1 : numel(v)-1
    if (v(j) ~= currentR && v(j) ~= v(j+1))
        if v(j) > 0
            NS{1} = [NS{1}, v(j)];
        else
            NS{2} = [NS{2}, v(j)];
        end
    end
end

if (numel(v) > 0 && v(end) ~= currentR )
    if v(end) > 0
        NS{1} = [NS{1}, v(end)];
    else
        NS{2} = [NS{2}, v(end)];
    end
end
numGM = numel(NS{2});
if numGM > 0
    tmpNS = zeros(numGM,5);
    pixelPos = zeros(numGM,2);
    cntrInROI = min(centroid, [1 1]+ND);
    for l = 1 : numGM
        [r,c] = find(currentBox == NS{2}(l));
        D = sqrt((r-cntrInROI(2)).^2 + (c-cntrInROI(1)).^2);
        [val, ind] =min(D);
        tmpNS(l,1) = - NS{2}(l); % germ mass label
        tmpNS(l,2) = val; % ditance of centriod and nearest pixel in boundary
        % first coordinate of nearest pixel or row
        tmpNS(l,3) = r(ind)+ max(1,centroid(2)-ND)- 1;
        % second coordinate of nearest pixel or column
        tmpNS(l,4) = c(ind)+max(1,centroid(1)-ND) - 1;
        pixelPos(l, :) = [c(ind) r(ind)];
        % weight sertoli-germMass connection by portion of set pixel in BW
        % per line length.
        [~, ~, cp] = improfile(BW, ...
            [centroid(1), tmpNS(l,4)], [centroid(2), tmpNS(l,3)]);
        setP = find(cp);
        tmpNS(l,5) = length(setP); %/length(cp);
    end
    % sort neighboring germ mass based on distance
    [~, sortTmpNS] = sort(tmpNS(:,2)); 
    NS{2} = tmpNS(sortTmpNS, : );
    pixelPos = pixelPos(sortTmpNS, :);
    NS{3} = checkBridgeRegion(currentR, cntrInROI, NS{2},pixelPos, currentBox);
end
end

function BR = checkBridgeRegion(currentR, centroid, GM, pixelPos, lbl)
%neighborhood distance of centroid that is applied to check direct conection
n= 5;
numGM = size(GM,1);
BR = cell(numGM,1);
[row col] = size(lbl);

for i = 1 : numGM
    % collecting some pixels in between regions
    connectingPixels = [];
    NearstP = pixelPos(i,:);
    [x, y, cp] = improfile(lbl, ...
        [centroid(1), NearstP(1)], [centroid(2), NearstP(2)]);
    connectingPixels = [connectingPixels; cp];
    cp = improfile(lbl, ...
        max([centroid(1)-n, NearstP(1)],1), [centroid(2), NearstP(2)]);
    connectingPixels = [connectingPixels; cp];
    cp = improfile(lbl, ...
        min([centroid(1)+n, NearstP(1)], col), [centroid(2), NearstP(2)]);
    connectingPixels = [connectingPixels; cp];
    cp = improfile(lbl, ...
        [centroid(1), NearstP(1)], max([centroid(2)-n, NearstP(2)], 1));
    connectingPixels = [connectingPixels; cp];
    cp = improfile(lbl, ...
        [centroid(1), NearstP(1)], min([centroid(2)+n, NearstP(2)], row));
    connectingPixels = sort([connectingPixels; cp]);
    
    for j = 1 : numel(connectingPixels)-1
        if (connectingPixels(j) ~= 0 && connectingPixels(j) ~= currentR ...
                && connectingPixels(j) ~= connectingPixels(j+1) ...
                && connectingPixels(j) ~= -GM(i,1) )
            BR{i} = [BR{i} connectingPixels(j)];
        end
    end
    if (connectingPixels(end) ~= 0 && connectingPixels(end) ~= currentR ...
            && connectingPixels(j) ~= -GM(i,1) )
        BR{i} = [BR{i} connectingPixels(end)];
    end
end
end


function directNS = directNeighbor(currentR, centroid, NS, centers, lbl, row, col)
%neighborhood distance of centroid that is applied to check direct conection
n= 5;

directNS = [];
for i = 1 : numel(NS)
    % collecting some pixels in between regions
    connectingPixels = [];
    NCentroid = centers(NS(i), :);
    [x, y, cp] = improfile(lbl, ...
        [centroid(1), NCentroid(1)], [centroid(2), NCentroid(2)]);
    connectingPixels = [connectingPixels; cp];
    cp = improfile(lbl, ...
        max([centroid(1)-n, NCentroid(1)-n],1), [centroid(2), NCentroid(2)]);
    connectingPixels = [connectingPixels; cp];
    cp = improfile(lbl, ...
        min([centroid(1)+n, NCentroid(1)+n], col), [centroid(2), NCentroid(2)]);
    connectingPixels = [connectingPixels; cp];
    cp = improfile(lbl, ...
        [centroid(1), NCentroid(1)], max([centroid(2)-n, NCentroid(2)-n], 1));
    connectingPixels = [connectingPixels; cp];
    cp = improfile(lbl, ...
        [centroid(1), NCentroid(1)], min([centroid(2)+n, NCentroid(2)+n], row));
    connectingPixels = [connectingPixels; cp];
    
    for j = 1 : numel(connectingPixels)
        if (connectingPixels(j) ~= 0 && connectingPixels(j) ~= currentR ...
                && connectingPixels(j) ~= NS(i))
            break;
        end
    end
    if j == numel(connectingPixels)
        directNS = [directNS, NS(i)];
    end
end
end


function np = separateNeighbor(nGermMass, NS, cntr, centers, closeGM)
np = repmat(struct('neg',{},'pos',{}),nGermMass,1);

% compute the line between Sertoli's center and Germ mass border and divide
% neighboring Sertoli to two classes negetive(left) and positive(right)
nNbSertoli = size(NS{1},2);
for i = 1 : nGermMass
    p = NS{2}(i,[4,3]);
    m= (cntr(1) - p(1))/(cntr(2) - p(2));
    b = cntr(1) - m*cntr(2);
    neg = [];
    negGMlbl = [];
    negGMDist = [];
    pos =[];
    posGMlbl = [];
    posGMDist = [];
    for j = 1 : nNbSertoli
        Ncntr = centers(NS{1}(j),:);
        if (Ncntr(1)- m *Ncntr(2))< b
            neg = [neg, j];
            [gmL gmD]=getGMlbl(NS{1}(j), closeGM);
            negGMlbl = [negGMlbl, gmL];
            negGMDist = [negGMDist, gmD];
        else
            pos= [pos, j];
            [gmL gmD]=getGMlbl(NS{1}(j), closeGM);
            posGMlbl = [posGMlbl, gmL];
            posGMDist = [posGMDist, gmD];
        end
    end
    np(i).neg = neg;
    np(i).negGMlbl = negGMlbl;
    np(i).negGMDist = negGMDist;
    np(i).pos = pos;
    np(i).posGMlbl = posGMlbl;
    np(i).posGMDist = posGMDist;

end
end

function [germlbl DD]= getGMlbl(sID,closeGM)
NS = closeGM{sID};
germlbl = 0;
DD = 0;
for j =1 : size((NS{2}),1)
    if isempty(NS{3}{j}) 
        germlbl = NS{2}(j,1);
        DD = NS{2}(j,2);
        return
    end
end
end