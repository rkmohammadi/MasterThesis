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

for f = 4:  -1 : 4
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
  %
    rawFeatures = [msrROI.Mean]',%; msrROI.MinVal; msrROI.Perimeter; msrROI.Size;...
%         msrROI.StdDev; msrROI.Skewness; msrROI.ExcessKurtosis]';
    nRawFeatures = normalizeFeature(rawFeatures,'Zero min');
    
    K = 4; %4
    opts = statset('Display','final');
    [idx,ctrs,sumdist] = kmeans(nRawFeatures,K,'dist','sqEuclidean', 'start', 'cluster',...
        'replicates',5,'EmptyAction','singleton','Options',opts);
   
    % choose regino with low gray value which showthe germ mass.
    germCluster = 1;
    if ctrs(2,1) < ctrs(1,1)
        germCluster = 2;
    end
    
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
        if c == germCluster
            title(['(' imgFile ') Cluster ' 'Germ Mass']);
        else
            title(['(' imgFile ') Cluster ' int2str(c) ' out of ' int2str(K)]);
        end
        set(distIm,'AlphaData',0.3);
    end
    
    % show region connectivity graph  
    germLbl = bwlabel(clusterdRgn == germCluster);
    stats = regionprops(germLbl, 'Centroid');
    neiborhood = 35;
    tic
    NRG = findNeighbors(germLbl, stats, neiborhood);
    toc
    % show graph
    figure, imshow(en_im), %lbl_bin_FRST),
    % title(['connected regions', fileName]);
    for i =1 : size(stats,1)
        NS = NRG{i};
        for j =1 : numel(NS)
            line([stats(i).Centroid(1), stats(NS(j)).Centroid(1)], ...
                [stats(i).Centroid(2), stats(NS(j)).Centroid(2)]);
        end
    end
end