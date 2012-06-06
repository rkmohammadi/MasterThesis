freqs = [0.5 1 2 3 4];
scales = [4 5];
ori = linspace(0,2,8)*pi/2;
ori = ori(1:end-1);

filterBank = FilterBank();
filterBank.Frequencies = freqs;
filterBank.Scales = scales;
filterBank.Orientations = ori;

filterBank.CreateFilterBank();

% filterBank.ShowFilters();
    
fileName =  '../../images/GATA-4/M3G4015.TIF';
Im = imread(fileName);
Im = imresize(Im,0.5,'lanczos2');
imPatch = rgb2gray(Im); %imread('lena.jpg'));

%imPatch = imPatch(1:100,1:100);

[filtersParams, responses] = filterBank.Convolve(double(imPatch));
%filterBank.ShowResponses(responses, filtersParams); 

%% cluster pixels
szR = size(responses);
szFiltOut = size(responses{1});

X =zeros(szFiltOut(1)*szFiltOut(2), szR);
for i=1 :szR
    temp = responses{i};
    X(:,i) = temp(:);
end

 K = 70;
opts = statset('Display','final');
[idx,ctrs,sumdist] = kmeans(X,K,'dist','sqEuclidean',...
'replicates',5,'Start','sample', 'EmptyAction','singleton','Options',opts);

%%
featureValues = zeros(1,length(freqs));
for i = 1:length(freqs);
    tmp = zeros(1,length(ori));
    for j = 1:length(ori)
        tmp(j) = mean(responses{(i-1)*length(ori)+j}(:));
    end
    featureValues(i) = mean(tmp);
end

% f = real(filterBank.FiltersValues{5});
% f = f - mean(f(:));
% f = f./max(abs(f(:)));
% f = uint8(f.*255+127);
% % imwrite(f,'GaborEx1.png');
%
% f = real(filterBank.FiltersValues{6});
% f = f - mean(f(:));
% f = f./max(abs(f(:)));
% f = uint8(f.*255+127);
% % imwrite(f,'GaborEx2.png');
%
% f = real(filterBank.FiltersValues{7});
% f = f - mean(f(:));
% f = f./max(abs(f(:)));
% f = uint8(f.*255+127);
% % imwrite(f,'GaborEx3.png');
%
% f = real(filterBank.FiltersValues{9});
% f = f - mean(f(:));
% f = f./max(abs(f(:)));
% f = uint8(f.*255+127);
% figure;imshow(f,[]);
% % imwrite(f,'GaborEx4.png');
