% manually define class of some region
% the file name are stored invariable images

images = {'../images/GATA-4/M3G4001.TIF'; ...
    '../images/GATA-4/M3G4004.TIF'; ...
    '../images/GATA-4/M3G4006.TIF'; ...
    '../images/GATA-4/M1G40019.TIF'; ...
    '../images/GATA-4/M1G40022.TIF'; ...
    '../images/GATA-4/M1G40036.TIF'};

szRange = [50 125;100 225;200 350;300 550];

for f = numel(images)-1 :  -1 : 1
    imgFile = images{f};
    [diplbl, en_im] = PreProcess(imgFile);
    sz = size(en_im);
    
    
    
    % apply measure to compute some properties for the region
    
    measurmentIDs = {'Size', 'Perimeter', 'Feret', 'Radius', ...%'ConvexPerimeter','Convexity',  ...
        'P2A', 'PodczeckShapes', 'Mean', 'StdDev', ...
        'Skewness', 'ExcessKurtosis', 'MaxVal', 'MinVal', 'Center', 'Inertia',...
        'Mu', 'MajorAxes'};
    lab = RGB2Lab(en_im);
    stretchedL = percentileStretch(lab(:,:,1), [0,1]);
    msr = measure(diplbl, mat2im(stretchedL),measurmentIDs, [],1);
    
    % sort msr based on Id
    msrId = msr.id;
    [~, sortId] = sort(msrId);
    %msr = msr(sortId);
    for rng = 1 : size(szRange,1)
        % array index of region in the current range
        cur = find((msr.Size > szRange(rng,1)) & (msr.Size < szRange(rng,2)));
        
        msrROI = zeros(size(cur));
        numR = numel(cur);
        % extract the msr of the current regions
        msrROI = msr(msrId(cur));
        
        %
        %     msrROI = msr; %msr([class1Region; class2Region]);
        rawFeatures = [msrROI.Mean; msrROI.MinVal; msrROI.Perimeter; ...%msrROI.Size;...
            msrROI.StdDev; msrROI.Skewness; msrROI.ExcessKurtosis]';
        % rawFeatures = [msrROI.Skewness; msrROI.ExcessKurtosis]';
        
        [fr, ~]=size(rawFeatures);
        nRawFeatures = (rawFeatures - ones(fr,1)*mean(rawFeatures)) ./ (ones(fr,1)*std(rawFeatures));
        
        K = 2; %4
        opts = statset('Display','final');
        [idx,ctrs,sumdist] = kmeans(nRawFeatures,K,'dist','sqEuclidean', 'start', 'cluster',...
            'replicates',5,'EmptyAction','singleton','Options',opts);
        
%         rId = msrROI.id;
%         sortedIdx = idx(sortId);
        double_lbl = double(diplbl);
        clusterdRgn= zeros(sz(1:2));
%         for row = 1 : sz(1)
%             for col = 1 : sz(2)
%                 if double_lbl(row,col)
%                     clusterdRgn(row,col) = sortedIdx(double_lbl(row,col));
%                 end
%             end
%         end
        %     for r =1 : length(rId)
        %         y(double_lbl == rId(r)) = idx(r);
        %     end
        
        for r =1 : numR
                 clusterdRgn(double_lbl == msrId(cur(r))) = idx(r);
        end
        
        sz = size(en_im);
        for c = 1: K
            figure,imshow(en_im);
            RGBDist = uint8(cat(3, zeros(sz(1:2)),(clusterdRgn == c)*255, ...
                zeros(sz(1:2))));
            
            hold on
            distIm = imshow(RGBDist);
            title(['(' imgFile ') Cluster ' int2str(c) ' out of ' int2str(K) ...
                ,'of size[' int2str(szRange(rng,1)) ',' int2str(szRange(rng,2)) ']']);
            set(distIm,'AlphaData',0.3);
        end
        
    end
    
end