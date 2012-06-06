function showSerotliAng(CloseGM, ang, centers, GMandSertoli, en_im, posNeg)
% function showSerotliAng(CloseGM, ang, centers, GMandSertoli, en_im, posNeg)
% Use this function to show connection line between srtoli around an
% germ-mass by considering the angle between the cell and one of its
% neighbors and the closest point on the germ-mass.
%  Input:
%   CloseGM, neighborhood information of the Sertoli cells.
%   ang, angle between Sertoli cells and the germ-mass.
%   centers, the centroids of the Sertolli cells.
%   GMandSertoli, the labeled image contains both germ-mass and Sertoli
%   cells.
%   en_im, the color image of the testes tissue.
%   posNeg, an structure to specify each neighbor of a Sertoli cell is in 
%   which side of its connection with the germ-mass.
% 
% See also findClosestGM, combineGMLumen, findNeighbors, makeAnglePositive.
% 
% Author: A.Rahim Kadkhodamohammadi (r.k.mohammadi@gmail.com)
% 01 June 2012 CBA, Uppsala University
%--------------------------------------------------------------

% internal arguments
cutOff = .5;%1.2;
minDist = 15;

%if nargin > 5
    [~, h] = overlayCluster(en_im, GMandSertoli, [], 1); 
%end
hold on


num = size(centers,1);
%%
if nargin == 5
    for i = 1 : num
        SNode = CloseGM{i};
        cntr = centers(i,:);
        neighborNum = size((SNode{1}),2);
        tmpAng = ang{i};
        %     tmpAng = and(tmpAng > 0, tmpAng<cutOff);
        % %     if i == 68
        % %         disp('here');
        % %     end
        if neighborNum && numel(tmpAng)
            tmpInd = or(tmpAng(:,1) < -cutOff, tmpAng(:,1)>cutOff);
            tmpAng(tmpInd,1) = 10000;
            [val, ind] =min(abs(tmpAng(:,1)));
            %         for j =1 : neighborNum
            if val < 1000 %find(tmpAng(j,:),1,'first')
                line([cntr(1), centers(SNode{1}(ind),1)], ...
                    [cntr(2), centers(SNode{1}(ind),2)]);
            end
            %         end
        end
    end
else
    for i = 1 : num
        SNode = CloseGM{i};
        cntr = centers(i,:);
        neighborNum = size((SNode{1}),2);
        tmpAng = ang{i};
        %     tmpAng = and(tmpAng > 0, tmpAng<cutOff);
        % %     if i == 68
        % %         disp('here');
        % %     end
        if neighborNum && numel(tmpAng)
            tmpInd = or(tmpAng(:,1) < -cutOff, tmpAng(:,1)>cutOff);
            tmpAng(tmpInd,1) = 10000;
            % check for negative side
            if SNode{2}(1,2) > minDist
                if ~isempty(posNeg{i}(1).neg)
                    [val, ind] =min(abs(tmpAng(posNeg{i}(1).neg,1)));
                    if val < 1000
                        rId = posNeg{i}(1).neg(ind);
                        line([cntr(1), centers(SNode{1}(rId),1)], ...
                            [cntr(2), centers(SNode{1}(rId),2)]);
                    end
                end
                if ~isempty(posNeg{i}(1).pos)
                    [val, ind] =min(abs(tmpAng(posNeg{i}(1).neg,1)));
                    if val < 1000
                        rId = posNeg{i}(1).neg(ind);
                        line([cntr(1), centers(SNode{1}(rId),1)], ...
                            [cntr(2), centers(SNode{1}(rId),2)]);
                    end
                end
            end
        end
    end
    
end

end