function normalFeatures = normalizeFeature(features, method)
% function normalFeatures = normalizeFeature(features, method)
% this function get a feature matrix that each row coressponde to an
% observation and normalize the features.
% If the number of input arguments is one normalize the feature to have
% zero means and unit standard deviation. 
% But if the second input is 'zero min', it normalize features to be
% positive.
%
% Author: A.Rahim Kadkhodamohammadi (r.k.mohammadi@gmail.com)
% 24/02/2012 CBA, Uppsala University
%--------------------------------------------------------------
if nargin < 1
    error('This function needs at least one input');
end

[fr, ~]=size(features);
if nargin == 1
    normalFeatures = (features - ones(fr,1)*mean(features)) ./ (ones(fr,1)*std(features));
elseif nargin ==2 && strcmpi(method,'zero min')
    normalFeatures = features - ones(fr,1)*min(features);
end

end