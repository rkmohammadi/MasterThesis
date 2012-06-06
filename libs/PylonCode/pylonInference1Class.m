function labels = pylonInference1Class(nBaseRegions, Tree, U, V, additionalUnaries)
%   finds the optimal configuration in the 1-class case
%   see the paper for more details.
%   V. Lempitsky, 2011
%   INPUT:
%       nBaseRegions -- number of leaves in the segmentation tree
%       Tree -- segmentation tree, that is (number non-leaf nodes)-by-2 list in the
%               format returned by matlab's 'cluster' function
%       U - 1x(number of all nodes/regions) arbitrary unary terms for
%               regions
%       V - 3x(number of edges between the leaf regions). Each column
%               should give the numbers of the two connected leaves and the
%               strength of the edge
%       additionalUnaries(optional) - additional unary terms for the x variables of leaf nodes only. Can
%               be used, e.g., to express Hamming loss during max-margin learning
%   OUTPUT:
%       labels -- output in 'x' notation (labels each node in the tree with
%               0,or 1). Use 'pylonConvertLabels' to obtain region
%               labels (i.e. set of regions that have been assigned to the foreground).

unary = zeros(2, numel(U));
unary(:,end) = [0; U(end)];
if nargin == 5 && ~isempty(additionalUnaries)
   unary(:,1:nBaseRegions) = additionalUnaries;
end

cluster

edgesBU = zeros(2, numel(U)-1);
pairwiseBU = zeros(4, numel(U)-1);

for i = 1:size(Tree,1)
  edgesBU(:,Tree(i,1)) = [Tree(i,1); i+nBaseRegions];
  pairwiseBU(3,Tree(i,1)) = U(Tree(i,1));
  edgesBU(:,Tree(i,2)) = [Tree(i,2); i+nBaseRegions];
  pairwiseBU(3,Tree(i,2)) = U(Tree(i,2));
end

pairwiseBU(2,:) = 1e10;

edgesBND = V(1:2,:);
pairwiseBND= [zeros(1,size(V,2)); V(3,:); V(3,:); zeros(1,size(V,2))];

labels = SolveQPBO(unary, [edgesBU edgesBND], [pairwiseBU pairwiseBND]);

