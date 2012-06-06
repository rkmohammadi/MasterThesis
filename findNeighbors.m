function NRG = findNeighbors(lbl, stats, d, img, thre)
% this function get the labeled image and distance for neighborhood and
% return neighboring region graph.
% if this function called with 5 argument, it will consider the gradient 
% infomation to decide weather connet or not. the function will consider
% this case if the input for stats is 2D array of region's centers. 
% 
% input:
%       lbl, labeled image. it could be double if stats is array of 
%           structure that contain field 'Centroid'. Otherwise it could 
%           be dip_image and stats should be an 2D array of center(2,n) 
%           that n is number of regions supposed to be gratter than 2 and 
%           each column coresspond to one region's center.
%       stats, -array of structure that at least has centriod of regions
%              - array of dip_measurement that at least has 'Center'
%                 (It could be computed by measure method of DIP_image). it
%                 is not recommended to call function with this input. it 
%                 will take a long time.   
%              - 2D array of size(2,n), it should contain the center of
%              region as it is computed by measure method of DIP_image.
%       d, neighborhood distance, an integer that specify the maximum 
%          distance between regions.
%       img, gray level image to compute the gradient of pixels in
%              between regions.
%       thre, threshold value
% output,
%       NRG : an array of structure that contain neighboring regions label
% 
% A.Rahim Kadkhodamohammadi (r.k.mohammadi@gmail.com)
% February 7 /2012 CBA, Uppsala University
%--------------------------------------------------------------------------

if (nargin < 3)
    error('this function needs at leat three arguments');
end

%turn off some warnings
wid = 'signal:findpeaks:largeMinPeakHeight';
warning('off', wid);

num = length(stats); % the number of regions
NRG = cell(num,1);
[row, col] = size(lbl);
switch class(stats)
    case 'struct'
        centers = round(cat(1, stats.Centroid));
        for i =1 : num
            centroid = centers(i,:);
            NS = neighboringSegment(i, centroid, row, col, lbl, d);
            if numel(NS) > 0
                directNS = directNeighbor(i, centroid, NS, centers, lbl, row, col);
                %  if you want to have pair of connected region uncomment following line
                %         for j= 1 : numel(directNS)
                %             NRG{directNS(j)} = [NRG{directNS(j)}, i];
                %         end
                NRG{i} = [NRG{i}, directNS];
            end
        end
    case 'dip_measurement'
        sId = stats.id;
        [~, sortId] = sort(sId);
        centers = stats.Center;
        centers = round(centers(:,sortId));
        for i =1 : num
            NS = neighboringSegmentDIP(i, centers(:,i), row, col, lbl, d);
            if numel(NS) > 0
                directNS = directNeighborDIP(i, centers(:,i), NS, centers, lbl, row, col);
                %  if you want to have pair of connected region uncomment following line
                %         for j= 1 : numel(directNS)
                %             NRG{directNS(j)} = [NRG{directNS(j)}, i];
                %         end
                NRG{i} = [NRG{i}, directNS];
            end
        end
    case 'double'
        if isa(lbl, 'dip_image')
            lbl = double(lbl);
            [row, col] = size(lbl);
        end
        if size(stats,1) < size(stats,2) % it is calculated by dipmeasure
            centers = round(stats' +1);
        end
        if nargin ==3
            for i =1 : num
                NS = neighboringSegment(i, centers(i,:), row, col, lbl, d);
                if numel(NS) > 0
                    directNS = directNeighbor(i, centers(i,:), NS, centers, lbl, row, col);
                    %  if you want to have pair of connected region uncomment following line
                    %         for j= 1 : numel(directNS)
                    %             NRG{directNS(j)} = [NRG{directNS(j)}, i];
                    %         end
                    NRG{i} = [NRG{i}, directNS];
                end
            end
        elseif nargin == 5
            for i =1 : num
                NS = neighboringSegment(i, centers(i,:), row, col, lbl, d);
                if numel(NS) > 0
                    directNS = directNeighborGrad(i, centers(i,:), NS, ...
                        centers, lbl, row, col, img, thre);
                    %  if you want to have pair of connected region uncomment following line
                    %         for j= 1 : numel(directNS)
                    %             NRG{directNS(j)} = [NRG{directNS(j)}, i];
                    %         end
                    NRG{i} = [NRG{i}, directNS];
                end
            end
        end
    otherwise
        error('stats should be one of these types: double, dip_measure or struct');
end
end

function NS = neighboringSegment(currentR, centroid, row, col, lbl, d)
[~, ~, v] = find(lbl(max(1,centroid(2)-d): min(centroid(2)+d, row), ...
    max(1,centroid(1)-d): min(centroid(1)+d, col))-currentR);
NS = [];
v = sort(v);
for j = 1 : numel(v)-1
    if (v(j) > 0 && v(j) ~= v(j+1))
        NS = [NS, v(j)+currentR];
    end
end

if (numel(v) > 0 && v(end)>0 ); NS=[NS, v(end)+currentR]; end
end

function NS = neighboringSegmentDIP(currentR, centroid, row, col, lbl, d)
[~, ~, v] = find(im2mat(lbl(max(0,centroid(1)-d): min(centroid(1)+d, row-1), ...
    max(0,centroid(2)-d): min(centroid(2)+d, col-1)))-currentR);
NS = [];
v = sort(v);
for j = 1 : numel(v)-1
    if (v(j) > 0 && v(j) ~= v(j+1))
        NS = [NS, v(j)+currentR];
    end
end

if (numel(v) > 0 && v(end)>0 ); NS=[NS, v(end)+currentR]; end
end

function directNS = directNeighborDIP(currentR, centroid, NS, centers, lbl, row, col)
%neighborhood distance of centroid that is applied to check direct conection
n= 5;

directNS = [];
for i = 1 : numel(NS)
    % collecting some pixels in between regions
    connectingPixels = [];
    NCentroid = centers(:,NS(i));
    [~, cp] = get_coords([centroid(1) centroid(2); ...
        NCentroid(1) NCentroid(2)],lbl);
    connectingPixels = [connectingPixels; cp];
    [~, cp] = get_coords([max(centroid(1)-n,0) centroid(2); ...
        max(NCentroid(1)-n, 0) NCentroid(2)],lbl);
    connectingPixels = [connectingPixels; cp];
    [~, cp] = get_coords([min(centroid(1)+n, row-1) centroid(2); ...
        min(NCentroid(1)+n, row-1) NCentroid(2)],lbl);
    connectingPixels = [connectingPixels; cp];
    [~, cp] = get_coords([centroid(1) max(centroid(2)-n,0); ...
        NCentroid(1) max(NCentroid(2)-n,0)],lbl);
    connectingPixels = [connectingPixels; cp];
    [~, cp] = get_coords([centroid(1) min(centroid(2)+n, col-1); ...
        NCentroid(1) min(NCentroid(2)+n, col-1)],lbl);
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

function directNS = directNeighborGrad(currentR, centroid, NS, centers, ...
    lbl, row, col, img, thre)
%neighborhood distance of centroid that is applied to check direct conection
n= 5;

directNS = [];
for i = 1 : numel(NS)
    % collecting some pixels in between regions
    connectingPixels = [];
    NCentroid = centers(NS(i), :);
    [x, y, c2cp] = improfile(lbl, ...
        [centroid(1), NCentroid(1)], [centroid(2), NCentroid(2)]);
    imgVal = impixel(img,x,y);
    imgVal = imgVal(:,1);
    connectingPixels = [connectingPixels; c2cp];
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
        if noEdge(currentR, NS(i), imgVal, c2cp, thre)
            directNS = [directNS, NS(i)];
        end
    end
end
end

function out = noEdge(curr, neighbor, imgVal, lblLine, thre)
% there are four scenarios for connecting line of two regions
% -first, the start from those regions,
% -second, the line start from one of them but not pass the other one
% -third, the line connect the region without crossing them
% -fourth, line start out regions but cross both or one of them

out = true;

g = gradient(imgVal);
nlbl = lblLine < 1;
if lblLine(1) == curr && lblLine(end) == neighbor
    f = find(nlbl, 1, 'first');
    e = find(nlbl, 1, 'last');
    if e-f > 5
        roi = abs(g(f+1:e-1));
        pks = findpeaks(roi, 'MINPEAKHEIGHT', thre);
        if numel(pks) > 0
            out =false;
        end
    end
elseif lblLine(1) == curr && lblLine(end) ~= neighbor
    f = find(nlbl, 1, 'first');
    e = find(lblLine == neighbor, 1, 'first');
    if isempty(e)
        e = length(g);
    end
    e = e-1;
    if e-f >3
        roi = abs(g(f-1:e));
        pks = findpeaks(roi, 'MINPEAKHEIGHT', thre);
        if numel(pks) > 0
            out =false;
        end
    end
elseif lblLine(1) ~= curr && lblLine(end) == neighbor
    f = find(lblLine == curr, 1, 'last');
    e = find(nlbl, 1, 'last');
    if isempty(f)
        f= 1;
    end
    f = f + 1;
    if e-f > 3
        roi = abs(g(f:e-1));
        pks = findpeaks(roi, 'MINPEAKHEIGHT', thre);
        if numel(pks) > 0
            out =false;
        end 
    end
else
    f = find(lblLine == curr, 1, 'last');
    e = find(lblLine == neighbor, 1, 'first');
    if isempty(f)
        f= 1;
    end
    f = f + 1;
    if isempty(e)
        e= length(g);
    end
    if e-f > 4
        roi = abs(g(f+1:e-1));
        pks = findpeaks(roi, 'MINPEAKHEIGHT', thre);
        if numel(pks) > 0
            out =false;
        end 
    end
end
end