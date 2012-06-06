% manually define class of some region
% the file name are stored invariable images
% this function can cluster the image region into region and then cluster
% each of them into two cluster, i.e it cluster the region into 2 and then
% cluster each region into 2 cluster again. As result we will have 4
% different cluster of regions.

images = {'../images/GATA-4/M3G4001.TIF'; ...
    '../images/GATA-4/M3G4004.TIF'; ...
    '../images/GATA-4/M3G4006.TIF'; ...
    '../images/GATA-4/M1G40019.TIF'; ...
    '../images/GATA-4/M1G40022.TIF'; ...
    '../images/GATA-4/M1G40036.TIF'};

for f = numel(images) :  -1 : 1
    imgFile = images{f};
    [diplbl, en_im] = PreProcess(imgFile);
    sz = size(en_im);
    
    
    
    % apply measure to compute some properties for the region
    
    measurmentIDs = {'Size', 'Perimeter', 'Feret', 'Radius', 'ConvexPerimeter', ...
        'P2A', 'PodczeckShapes', 'Convexity', 'Mean', 'StdDev', ...
        'Skewness', 'ExcessKurtosis', 'MaxVal', 'MinVal', 'Center', 'Inertia',...
        'Mu', 'MajorAxes'};
    lab = RGB2Lab(en_im);
    stretchedL = percentileStretch(lab(:,:,1), [0,1]);
    msr = measure(diplbl, mat2im(stretchedL),measurmentIDs, [],1);
    
    % sort msr based on Id
    tmp = msr.id;
    [~, sortId] = sort(tmp);
    %msr = msr(sortId);
    
    %
    msrROI = msr; %msr([class1Region; class2Region]);
%     rawFeatures = [msrROI.Mean; msrROI.MinVal; msrROI.Perimeter; msrROI.Size;...
%         msrROI.StdDev; msrROI.Skewness; msrROI.ExcessKurtosis]';
    % rawFeatures = [msrROI.Skewness; msrROI.ExcessKurtosis]';
  %%
    rawFeatures = [msrROI.Mean; msrROI.MinVal; msrROI.Perimeter; msrROI.Size;...
        msrROI.StdDev; msrROI.Skewness; msrROI.ExcessKurtosis]';
    nRawFeatures = normalizeFeature(rawFeatures);
    
    K = 2; %4
    opts = statset('Display','final');
    [idx,ctrs,sumdist] = kmeans(nRawFeatures,K,'dist','sqEuclidean', 'start', 'cluster',...
        'replicates',5,'EmptyAction','singleton','Options',opts);
    
    rId = msrROI.id;
    sortedIdx = idx(sortId);
    double_lbl = double(diplbl);
    clusterdRgn= zeros(sz(1:2));
    for row = 1 : sz(1)
        for col = 1 : sz(2)
            if double_lbl(row,col)
                clusterdRgn(row,col) = sortedIdx(double_lbl(row,col));
            end
        end
    end
%     for r =1 : length(rId)
%         y(double_lbl == rId(r)) = idx(r);
%     end
    
    sz = size(en_im);
    for c = 1: K
        figure,imshow(en_im);
        RGBDist = uint8(cat(3, zeros(sz(1:2)),(clusterdRgn == c)*255, ...
            zeros(sz(1:2))));
        
        hold on
        distIm = imshow(RGBDist); 
        title(['(' imgFile ') Cluster ' int2str(c) ' out of ' int2str(K)]);
        set(distIm,'AlphaData',0.3);
    end
    %%
    % cluster the region of each cluster into two cluster based on their
    % intensity features
    c1= find(sortedIdx == 1); %region indices of the first cluster
    c2 = find(sortedIdx == 2); %region indices of the second cluster
    msrC1 = msr(c1);
    msrC2 = msr(c2);
    
    c1Features = [msrC1.Mean; msrC1.MinVal; %msrC1.Perimeter; msrC1.Size;...
        msrC1.StdDev; msrC1.Skewness; msrC1.ExcessKurtosis]';
    c2Features = [msrC2.Mean; msrC2.MinVal; %msrC2.Perimeter; msrC2.Size;...
        msrC2.StdDev; msrC2.Skewness; msrC2.ExcessKurtosis]';
    
    nC1Features = normalizeFeature(c1Features);
    nC2Features = normalizeFeature(c2Features);
    
    for subC = 1 : 2
        curFeature = eval(['nC' num2str(subC) 'Features']);
        opts = statset('Display','final');
        [idx,ctrs,sumdist] = kmeans(curFeature,K,'dist','sqEuclidean', 'start', 'cluster',...
            'replicates',5,'EmptyAction','singleton','Options',opts);
        
        clusterdRgn= zeros(sz(1:2));
        subId = eval(['c' num2str(subC)]);
        for r =1 : length(subId)
            clusterdRgn(double_lbl == subId(r)) = idx(r);
        end
        
        sz = size(en_im);
        for c = 1: K
            figure,imshow(en_im);
            RGBDist = uint8(cat(3, zeros(sz(1:2)),(clusterdRgn == c)*255, ...
                zeros(sz(1:2))));
            
            hold on
            distIm = imshow(RGBDist);
            title(['(' imgFile ') Cluster ' int2str(c) ' out of ' int2str(K)]);
            set(distIm,'AlphaData',0.3);
        end
    end
end