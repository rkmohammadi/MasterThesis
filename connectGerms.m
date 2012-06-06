function [filledRegion, filledRegion2] = connectGerms(BW, Sertoli)
% function [filledRegion, filledRegion2] = connectGerms(BW, Sertoli)
% It connect all germ cells and fill small background holes between the
% cells connections to get a concrete regions as the germ-mass regions.
% If this function called with only one input argument it will only
% implement the method that is described in section 4.1.3 "Connect germ
% cells". But if we provide the second input it will remove connections
% that pass through Sertil cells.
% Input:
%   BW, a binary image contains the germ cells.
%   Sertoli (obtional), a binary image contains the Sertoli cells.
% Output:
%   filledRegion, a binary image shows the germ-mass regions that is
%   constructed excatly as the described method in the report. 
%   filledRegion2, The germ-mass regions that are constructed with
%   considering Sertoli cells.
% 
% Author: A.Rahim Kadkhodamohammadi (r.k.mohammadi@gmail.com)
% 01 June 2012 CBA, Uppsala University
%--------------------------------------------------------------

if nargout > nargin
    error('Please check the function for correct number of input and ouput');
end

% internal parameter 
MaxDist = 55; % maximum neighborhood distance 
MaxSizeofHoles = 2000;
SEsz2removeConn = 15; % the size of the stuctural element to remove 
% inter-tubules connections 

lbl = bwlabel(BW,4);
stats = regionprops(lbl, 'Centroid', 'PixelList');
NRG = findNeighbors(lbl, stats, MaxDist);
connectedNeighbors = connectNeighbors(NRG, stats, BW); 
[lblConnected num] = bwlabel(connectedNeighbors,8);

filledRegion = zeros(size(BW));
for r = 1 : num
    R = (lblConnected == r);
%     R = CheckImageBoundary(R);
    fillR = imfill(R,'holes');
    filledRegion = or(fillR, filledRegion);
end

if nargin == 2 && nargout == 2
    [lblConnected num] = bwlabel(logical(and(connectedNeighbors, not(Sertoli))),8);
    filledRegion2 = zeros(size(BW));
    for r = 1 : num
        R = (lblConnected == r);
        fillR = imfill(R,'holes');
        filledRegion2 = or(fillR, filledRegion2);
    end
end

% remove one pixel lines that connect different regions.
% % % % % % % filledRegion = erosion(filledRegion, 3, 'elliptic');
% apply tophat transform to remove small regions that connect big regions
% that are detected as germ mass.
% % % % % % topF = tophat(filledRegion, 15);
% % % % % % filledRegion = and(filledRegion,not(topF>0));
holes = filledRegion & ~connectedNeighbors;
holes = imerode(holes,strel('disk',1));
bigHoles = bwareaopen(holes, MaxSizeofHoles, 4);
filledRegion = opening(filledRegion & ~bigHoles, SEsz2removeConn);
filledRegion = imfill(im2mat(filledRegion), 'holes');
end


function BW = connectNeighbors(NRG, stats, BW)

for r = 1 : length(NRG)
    c = round(stats(r).Centroid);
    if BW(c(2),c(1)) == 0
        c = stats(r).PixelList(end,:);   
    end
    NR = NRG{r};
    for i = 1 : numel(NR)
        cn = round(stats(NR(i)).Centroid);
        if BW(cn(2),cn(1)) == 0
            cn = stats(NR(i)).PixelList(1,:);
        end
        [x, y, ~] = improfile(BW, [c(1), cn(1)], [c(2), cn(2)]);
        x = round(x);
        y = round(y);
        for p = 1 : length(x)
            BW(y(p),x(p)) = 1;
        end
    end
end
end

function R = CheckImageBoundary(R)
sz = size(R);
% put boundary of image in Boundary
Boundary = int8([R(:,1)' R(sz(1), :) , R(sz(1):-1:1,sz(2))' R(1, sz(2):-1:1)]);
Boundary(2:end) = Boundary(2:end)-Boundary(1:end-1);
Boundary(1) = 0;
[~, ind, v] = find(Boundary);
if ~isempty(ind) && length(ind) > 3
   if v(1) == 1
       CP = findConnectingPoints(ind, v, numel(Boundary));
   else % first point is -1 then the begining of region is at the end of Boundary
       tmpV = v(end);
       tmpInd = ind(end);
       v(2:end) = v(1:end-1);
       v(1) = tmpV;
       ind(2:end) = ind(1: end-1);
       ind(1) = tmpInd;
       CP = findConnectingPoints(ind, v, numel(Boundary));
   end
   % change R for pixel indicies that are in CP. 
   newBoundary = zeros(size(Boundary));
   start = 1;
   if CP(1,1) > CP(2,1)
       newBoundary(1:CP(1,1)) = 1;
       newBoundary(end:-1:CP(2,1)) = 1;
       start = 2;
   end
   for ps = start : size(CP,2)
       newBoundary(CP(1,ps):CP(2,ps)) = 1;
   end
   % copy new value into boundary of image
   R(:,1) = newBoundary(1:sz(1))';
   R(sz(1),:) = newBoundary((sz(1)+1):(sz(1)+sz(2)));
   R(sz(1):-1:1,sz(2)) = newBoundary((sz(1)+sz(2)+1):(2*sz(1)+sz(2)));
   R(1, sz(2):-1:1) = newBoundary((2*sz(1)+sz(2)+1):end);
end

end

function CP = findConnectingPoints(ind, v, numOfPixel)
num = length(ind);
% if mod(num,4)
%     error('The number of region points that touch image boundary must be even!!');
% end

CP = zeros(2,num/4);
j = 1;
i = 1;
while i < num
    if v(i+1) ~= -1 || v(i+2) ~= 1
        error('The squence of point should be as 1,-1 ...');
    else
        d1 = ind(i+2)-ind(i+1); % forward distance
    end
    if i == 1, p = num; else p =i-1; end
    if v(p) ~= -1 || v(i) ~= 1
        error('The squence of point should be as 1,-1 ...');
    else
        if i == 1
            d2 = numOfPixel-ind(p) + ind(i);
        else
            d2 =ind(i)-ind(p);% backward distance
        end
    end
    
    if d1 < d2
        i = 1 + 4;
        CP(1,j) = ind(p) - 1;
        CP(2,j) = ind(i);
    else
        i = i + 6;
        CP(1,j) = ind(i+1) - 1;
        CP(2,j) = ind(i+2);
    end
    j = j + 1;
end
end