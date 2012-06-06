function labels = pylonInferenceMultiClass(nBaseRegions, Tree, U, V, additionalUnaries)
%   finds the optimal configuration in the multi-class case for M labels)
%   see the paper for more details.
%   V. Lempitsky, 2011
%   INPUT:
%       nBaseRegions -- number of leaves in the segmentation tree
%       Tree -- segmentation tree, that is (number non-leaf nodes)-by-2 list in the
%               format returned by matlab's 'cluster' function
%       U - Mx(number of all nodes/regions) arbitrary unary terms for
%               regions
%       V - 3x(number of edges between the leaf regions). Each column
%               should give the numbers of the two connected leaves and the
%               strength of the edge
%       additionalUnaries(optional) - additional unary terms for the x variables of leaf nodes only. Can
%               be used e.g. to express Hamming loss during max-margin learning
%   OUTPUT:
%       labels -- output in 'x' notation (labels each node in the tree with
%               0,1..M). Use 'pylonConvertLabels' to obtain region
%               labels (i.e. set of regions that have been assigned to each class).

nLabels = size(U,1);
nRegions = size(U,2);
[dummy, t] = min(U(:,1:nBaseRegions));
labels = [t zeros(1,nRegions-nBaseRegions)];
infty = 1e10;

U_ = [infty*ones(1,nRegions); U];

clip = find(V(3,:) < 0);
V(3,clip) = 0;

for iter = 1:100
    oldLabeling = labels;
    for alpha = 1:nLabels
        Ualpha = [U_(sub2ind(size(U_),labels+1,1:nRegions)); U(alpha,:)];
        Valpha = [V(1:2,:);...
            labels(V(1,:)) ~= labels(V(2,:)); labels(V(1,:)) ~= alpha; ...
            labels(V(2,:)) ~= alpha; zeros(1,size(V,2))];
        Valpha(3:6,:) = bsxfun(@times,Valpha(3:6,:),V(3,:));
        if nargin <= 4 || isempty(additionalUnaries)
            Lalpha = [];
        else
            Lalpha = [additionalUnaries(sub2ind(size(U),...
                labels(1:nBaseRegions),1:nBaseRegions)); additionalUnaries(alpha,:)];
        end      
        
        alphaLab = pylonInference2Class(nBaseRegions,Tree,Ualpha,Valpha, Lalpha);
        labels(alphaLab == 0) = 0;
        labels(alphaLab == 2) = alpha;       
    end
    if all(oldLabeling == labels)
        break
    end
end
