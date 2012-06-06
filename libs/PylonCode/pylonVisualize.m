function outImage = pylonVisualize(xlabels, nBaseRegions, baseRegionMap, Tree, palette)
%   finds the optimal configuration in the multi-class case for M labels)
%   see the paper for more details.
%   V. Lempitsky, 2011
%   INPUT:
%       xlabels -- labels in 'x' format
%       baseRegionMap -- a map of the same size as the input image, where
%                each pixel is assigned the number of the leaf region
%       nBaseRegions -- number of leaves in the segmentation tree
%       Tree -- segmentation tree, that is (number non-leaf nodes)-by-2 list in the
%               format returned by matlab's 'cluster' function
%       palette -- each row should give an RGB color for each label ([1.0 0.0 0.0] = red)
%   OUTPUT:
%       labels -- output in 'x' notation (labels each node in the tree with
%               0,1..M). Use 'pylonConvertLabels' to obtain region
%               labels (i.e. set of regions that have been assigned to each class).
    
pixelLabels = xlabels(baseRegionMap);

Rout = reshape(palette(pixelLabels,1),size(pixelLabels))*255;
Gout = reshape(palette(pixelLabels,2),size(pixelLabels))*255;
Bout = reshape(palette(pixelLabels,3),size(pixelLabels))*255;

%the tricky bit is to visualize region boundaries
%
regionLabels = -ones(size(Tree,1)+nBaseRegions,1);
for i = size(Tree,1):-1:1
  if xlabels(i+nBaseRegions) == 0
      continue;
  end
  if regionLabels(i+nBaseRegions) == -1
     regionLabels(i+nBaseRegions) = i+nBaseRegions;
  end
  regionLabels(Tree(i,1)) = regionLabels(i+nBaseRegions);
  regionLabels(Tree(i,2)) = regionLabels(i+nBaseRegions);
end
for i = 1:nBaseRegions
  if regionLabels(i) == -1
      regionLabels(i) = i;
  end
end    
Rout = Rout*0.7;
Gout = Gout*0.7;
Bout = Bout*0.7;
regionMask = regionLabels(baseRegionMap);

pL = unique(regionLabels(:));    
for i=1:numel(pL)
    if pL(i) == -1
        continue;
    end
    t = regionMask == pL(i);
    t = t-imerode(t,strel('diamond',1));
    Rout(t > 0) = 255;
    Gout(t > 0) = 255;
    Bout(t > 0) = 255;
end    

outImage = uint8(cat(3,Rout,Gout,Bout));
