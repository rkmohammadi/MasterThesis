function labelsOut = pylonConvertLabels(labelsIn, Tree, nBaseRegions)
%   converts labels from "x-representation" to "f-representation"
%   see the paper for more details
%   V. Lempitsky, 2011

labelsOut = zeros(2*size(Tree,1)+1,1);

if labelsIn(end) ~= 0
    labelsOut(end) = labelsIn(end);
end

for i = 1:size(Tree,1)
  if labelsIn(i+nBaseRegions) == 0 && labelsIn(Tree(i,1)) ~= 0
      labelsOut(Tree(i,1)) = labelsIn(Tree(i,1));
  end
  if labelsIn(i+nBaseRegions) == 0 && labelsIn(Tree(i,2)) ~= 0
      labelsOut(Tree(i,2)) = labelsIn(Tree(i,2));
  end 
end

