function out = makeAnglePositive(ang, closeGM)
% function out = makeAnglePositive(ang, closeGM)
% this function get two cell array contains the ang between the Sertoli
% cells and closeGM that represent the neighborhood infomation of Sertoli
% cells.
% 
% Author: A.Rahim Kadkhodamohammadi (r.k.mohammadi@gmail.com)
% 01 June 2012 CBA, Uppsala University
%--------------------------------------------------------------
num = numel(ang);
out = cell(num,1);
for i = 1 : num
    Node = closeGM{i};
    if ~isempty(ang{i})
        A = ang{i};
        for s = 1 : size(Node{1},2)
            sId = Node{1}(s);
            A2 = ang{sId};
            idofCurrent = find(closeGM{sId}{1}==i); % id of current Sertoli in another one
            if ~isempty(closeGM{sId}{2})
                for g = 1 : size(Node{2},1)
                    gId = find(closeGM{sId}{2}(:,1)==Node{2}(g,1));
                    if ~isempty(gId)
                        A(s,g) = mean([abs(A(s,g)), abs(A2(idofCurrent, gId))]);
                    else
                        A(s,g) = abs(A(s,g));
                    end
                end
            end
        end
        out{i} = A;
        
    end
end
end