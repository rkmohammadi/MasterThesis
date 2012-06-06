function labels = pylonInference2Class(nBaseRegions, Tree, U, V, additionalUnaries)
%   finds the optimal configuration in the 2-class case
%   see the paper for more details.
%   V. Lempitsky, 2011
%   INPUT:
%       nBaseRegions -- number of leaves in the segmentation tree
%       Tree -- segmentation tree, that is (number non-leaf nodes)-by-2 list in the
%               format returned by matlab's 'cluster' function
%       U - 2x(number of all nodes/regions) arbitrary unary terms for
%               regions (for the labels 1 and 2)
%       V - either 3x(number of edges between the leaf regions). Each column
%               should give the numbers of the two connected leaves and the
%               strength of the edge
%           or 6x(number of edges between the leaf regions). Each column
%               should give the numbers of the two connected leaves and 
%               pairwise potential [E00; E01; E10; E11] (which should be
%               submodular, otherwise optimality is not guaranteed)
%       additionalUnaries(optional) - additional unary terms for the x variables of leaf nodes only. Can
%               be used, e.g., to express Hamming loss during max-margin learning
%   OUTPUT:
%       labels -- output in 'x' notation (labels each node in the tree with
%               0,1,or 2). Use 'pylonConvertLabels' to obtain region labels.

unary = cell(1,2);
edgesBU = cell(1,2);
pairwiseBU = cell(1,2);
edgesBND = cell(1,2);
pairwiseBND = cell(1,2);

nRegions = size(U,2);

infty = 1e10;

for varset = 1:2
    
    unary{varset} = zeros(2, size(U,2));
    unary{varset}(:,end) = [0; U(varset,end)];
    if nargin == 5 && ~isempty(additionalUnaries)
       unary{varset}(2,1:nBaseRegions) = additionalUnaries(varset,:);
    end
        

    edgesBU{varset} = zeros(2, size(U,2)-1);
    pairwiseBU{varset} = zeros(4, size(U,2)-1);

    for i = 1:size(Tree,1)
      edgesBU{varset}(:,Tree(i,1)) = [Tree(i,1); i+nBaseRegions]+nRegions*(varset-1);
      pairwiseBU{varset}(3,Tree(i,1)) = U(varset,Tree(i,1));
      edgesBU{varset}(:,Tree(i,2)) = [Tree(i,2); i+nBaseRegions]+nRegions*(varset-1);
      pairwiseBU{varset}(3,Tree(i,2)) = U(varset,Tree(i,2));
    end

    pairwiseBU{varset}(2,:) = infty;

    if varset == 2      
        edgesBND{varset} = V(1:2,:)+nRegions*(varset-1);
        if size(V,1) == 3
            pairwiseBND{varset} = [zeros(1,size(V,2)); V(3,:); V(3,:); zeros(1,size(V,2))];
        else
            pairwiseBND{varset} = [V(3:6,:)];
        end
    end
end

edgesConsist = [1:nRegions; nRegions+1:2*nRegions];
pairwiseConsist = [zeros(3,nRegions); ones(1,nRegions)*infty];
pairwiseConsist(1, 1:nBaseRegions) = infty;

out = SolveQPBO([unary{1} unary{2}], [edgesBU{1} edgesBU{2} edgesBND{2} edgesConsist],...
    [pairwiseBU{1} pairwiseBU{2} pairwiseBND{2} pairwiseConsist]);

if any(out < 0)
    warning('QPBO returned unlabeled nodes!');
end

labels = zeros(nRegions,1);
labels( out(1:nRegions) == 1) = 1;
labels( out(nRegions+1:end) == 1) = 2;