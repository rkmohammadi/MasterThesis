function ang = regionAngle(SMG, centers)
% function ang = regionAngle(SMG, centers)
%  this function get the Sertoli cells and the neighborhood information and
%  their centroid to compute ang between Sertoli cells and germ-masses.
%  Input: 
%   SMG, the neighborhood infomation of the Sertoli cells that can be
%   colected using function findCloseGM.
%   centers, a 2D array that eahc rows represent the centroids of the cells. 
%  Ouoput:
%   ang, a cell array contains the angle between the cell and each of its
%   neighbors according to the closest point in germ-mass. 

% 
% Author: A.Rahim Kadkhodamohammadi (r.k.mohammadi@gmail.com)
% 01 June 2012 CBA, Uppsala University
%--------------------------------------------------------------

init = 10000;
num = numel(SMG);
ang = cell(num,1);
for i = 1 : num
    cntr = centers(i,:);
    gNode = SMG{i};
    numSer = size(gNode{1},2);
    numGerm = size(gNode{2},1);
    if numSer && numGerm 
        a = ones(numSer,numGerm)*init;
        for s = 1 : numSer
            neighborSS = SMG{gNode{1}(s)}{2};
            cntrN = centers( gNode{1}(s), :);
            if ~isempty(neighborSS)
                for g =1 : numGerm
                    ind = find(neighborSS(:,1) == gNode{2}(g,1));%,1,'first');
                    if ~isempty(ind)
                        %reverse the sotred indicies
                        a(s,g) = OcallaghanAng(cntr,cntrN,neighborSS(ind,[4,3]));
                    end
                end
            end
        end
        ang{i} = a;
    end
        
end
end


function d = OcallaghanAng(cntr,cntrN, p)
x= norm(cntrN - p);
y= norm(cntr - cntrN);
z= norm(cntr - p);
d= (x^2 + y^2 - z^2)/(x*y);
end