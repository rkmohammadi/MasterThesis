function [UnaryT, pairwiseC] = computeDataCost(in1, in2, method)
% function [UnaryT, pairwiseC] = computeDataCost(in1, in2, method)
% This function set up costs for graph cuts. It can use two different
% method to compute data cost. The input arguments are changed based on the
% method that used to compute costs.
% First, Distance2clustercenter(default method):
%   Distance to cluster center, it will compute the distance of each cell
%   to the centroid of cluster in the feature space. This centroid can be
%   calculated by K-means clustering algorithm. 
%   Note that in this case the function can only produce one output cost
%   for data cost(UnaryT).
%   In this case the first input is the centroids of the clusters and the
%   second one is 2D matrix that eahc row of the matrix represents a cell
%   in the feature space.
%   The output is 2D matrix which each row correspond to a cell and each
%   column has the inverse distance of the cells' features to the centroid
%   of the cluster.
%  Second, RatioOfDifferentCluster:
%   Ratio of different Cluster, the spacial relashenship of the cells are
%   considered to compute both data costs and smoothness costs, which are
%   retuned by UneryT and pairwiseC, respectively.
%   The uneryT is caluted based on the density of clusters in the
%   nieghborhood of the cells and pairwiseC is computed based on the
%   dnesity of the cluster in the neighborhood of the cluster and the
%   cluster label of the cells. More information can be found in section
%   4.1.2 Clustering Cells in my report.
%   The inputs are:
%       in1, a cells arrey that each element of the contains the index of
%       all neighboring cells of the current cell;
%       in2, a 1D array that eahc element contains the cluster id of its
%       corresponding cell;
%       method, the is a string and should be "RatioOfDifferentCluster".
%   The ouptuts are:
%       UnaryT, a 2D array that each row is data cost for cell;
%       pairwiseC, 2D array, the size of eahc dimension of the array is 
%       equal to the number of cells. Zeros means no connection between the
%       cells. High value for each entry show high connectivity between the
%       cells.
% 
% Author: A.Rahim Kadkhodamohammadi (r.k.mohammadi@gmail.com)
% 01 June 2012 CBA, Uppsala University
%--------------------------------------------------------------

if nargin == 2
    method = 'Distance2clustercenter';
end
if strcmpi(method, 'Distance2clustercenter')
    [numLbls c1] = size(in1); % cluster center 
    [numRegions c2] = size(in2);% feature matrix of regions
    if numLbls < 2
        error('You should at least have two regions!!!');
    end
    if c1 ~= c2
        error('The number of column in first and second input should agree!!!');
    end
    UnaryT = zeros(numLbls, numRegions);
    for l = 1 : numLbls
        for r = 1: numRegions
            UnaryT(l,r) = 1./norm(in1(l,:) - in2(r,:) ); % changed to inverse first ****
        end
    end
    if nargout == 2
        error('Can not compute pairwise cost with the provided input');
    end
elseif strcmpi(method, 'RatioOfDifferentCluster')
    % in1 NRG graph of neighborhood region
    % in2 cluster id of each region
    numLbls = max(in2(:));
    numRegions = numel(in1);
    UnaryT = zeros(numLbls, numRegions);
    for r = 1 : numRegions
        NR = in1{r};
        for i = 1: length(NR)
            UnaryT(in2(NR(i)),r) = UnaryT(in2(NR(i)),r) +1;
            UnaryT(in2(r),NR(i)) = UnaryT(in2(r),NR(i)) +1;
        end
    end
    sumUnaryT = sum(UnaryT); 
    s = (max(sumUnaryT)) ./ sumUnaryT;% to have the number of label in [0 20]
    for l =1 : numLbls
        UnaryT(l,:) = UnaryT(l,:) ./ s;
    end
    if nargout == 2
        pairwiseC = weightEdge(UnaryT,in1, in2, numRegions);
    end
end
end

function pairwiseC = weightEdge(UT, NRG, clusterId, numRegions)

pairwiseC = zeros(numRegions, numRegions);
for r = 1 : numRegions
    NR = NRG{r};
    for i = 1 : numel(NR)
        % if the label of the region is the same then max of unary of this
        % label in the regions otherwise minimum of strech of eahc region
        % to other one
        if clusterId(r) == clusterId(NR(i))
            pairwiseC(r, NR(i)) = 2*max(UT(clusterId(NR(i)),r), UT(clusterId(r), NR(i)));
        else
            pairwiseC(r, NR(i)) = 2*max(UT(clusterId(NR(i)),r), UT(clusterId(r), NR(i))); % it minimum first****
        end
%         pairwiseC(NR(i), r) = UT(clusterId(r), NR(i)); %set neighborhood for both side
    end
end

end