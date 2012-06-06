function w = weightR(NRG, features)
% function w = weightR(NRG, features)
% calculate the distance between neiboring region 
% input:
%       NRG, cell that contain the neiboring region indices 
%       features, feaeture matrix that contain feature for each region. the
%       size of feature is n by m where n is the number of region and equal
%       to the number of element of 'NRG', m is the number of features.
% output:
%       w, cell contain distance between each region and its neighboring
%       regions. each element of w is an array that its row corissponds to
%       neighboring region and column corissponds to features`
% 
% See also findClosestGM.
% A.Rahim Kadkhodamohammadi (r.k.mohammadi@gmail.com)
% March 6 /2012, CBA Uppsala University
%--------------------------------------------------------------------------


% the number of regions and features for each region
[numR, numF] = size(features);
w = cell(numR, 1);
for r = 1 : numR
    NR = NRG{r};
    numNR = numel(NR);
    currentRF = features(r,:);
    currentW = zeros(numNR, numF);
    
    for i = 1 : numNR
        currentW(i, :) = abs(currentRF - features(NR(i), :));
    end
    
    w{r} = currentW;
end

end