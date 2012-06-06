clear;

filterBank = FilterBank();
filterBank.CreateFilterBank();
% or:filterBank.CreateFilterBank(FilterBank.CREATE_CROSSPRODUCT);

% or:filterBank.CreateFilterBank(FilterBank.CREATE_SEQUENTIAL); %erase some
% values in parameters (scale, freq... -> FilterBank)to produce the same
% number of values per parameter (otherwise it will produce error)

%******************Figure 1: Filters - realParts; imshow (default) ********************
filterBank.ShowFilters();

%******************Figure 2: Filters - realParts; plot (user defined) ********************
featExtractor = @(x) GaborKernel.GetImagParts(x);
displayFunction = @(x) surf(x);
filterBank.ShowFilters(featExtractor, displayFunction);
% just put '[]' for default -> eg. filterBank.ShowFilters([], displayFunction);


im = imread('lena.jpg'); im = rgb2gray(im); im = double(im);

%******************Figure 3: Responses - amplitudes (default) ********************
[filtersParams, responses] = filterBank.Convolve(im);
FilterBank.ShowResponses(responses, filtersParams);

%******************Figure 4: Responses - phases (user defined) ********************
featureExtractor = @(x) GaborKernel.GetPhases(x);
[filtersParams, responses] = filterBank.Convolve(im, featureExtractor);
FilterBank.ShowResponses(responses, filtersParams);