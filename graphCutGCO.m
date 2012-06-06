images = {'../images/GATA-4/M3G4001.TIF'; ...
    '../images/GATA-4/M3G4004.TIF'; ...
    '../images/GATA-4/M3G4006.TIF'; ...
    '../images/GATA-4/M1G40019.TIF'; ...
    '../images/GATA-4/M1G40022.TIF'; ...
    '../images/GATA-4/M1G40036.TIF'};



% imgFile = '../images/GATA-4/M1G40019.TIF';
folder = '../images/GATA-4/';
files = dir([ folder '*.TIF']);
for fi = 2 : 2% numel(files)
% imgFile = '../images/GATA-4/M1G40022.TIF';%images{fi};
imgFile = [folder files(fi).name];
%imgFile = '../images/GATA-4/M3G4015.TIF';
outFile = ['../../tmpResult' imgFile(17:end-3) 'tif'];
[diplbl, en_im tt, edg] = PreProcess(imgFile);
sz = size(en_im);

% apply measure to compute some properties for the region

measurmentIDs = {'Size', 'Perimeter', 'Mean', 'StdDev', ...
    'Skewness', 'MaxVal', 'MinVal', 'Center', 'Inertia'};
% {'Size', 'Perimeter', 'Feret', 'Radius', 'ConvexPerimeter', ...
%     'P2A', 'PodczeckShapes', 'Convexity', 'Mean', 'StdDev', ...
%     'Skewness', 'ExcessKurtosis', 'MaxVal', 'MinVal', 'Center', 'Inertia',...
%     'Mu', 'MajorAxes'};
lab = RGB2Lab(en_im);
%%
stretchedL = percentileStretch(lab(:,:,1), [0,1]);
msr = measure(diplbl, mat2im(stretchedL),measurmentIDs, [],1);

% sort msr based on Id
tmp = msr.id;
[~, sortId] = sort(tmp);
%msr = msr(sortId);

%
msrROI = msr;
rawFeatures = [msrROI.Mean; msrROI.MinVal; msrROI.StdDev; msrROI.Skewness]';% msrROI.Perimeter; msrROI.Size;...
        %msrROI.StdDev; msrROI.ExcessKurtosis]';
rawFeatures = rawFeatures(sortId,:);
nRawFeatures = normalizeFeature(rawFeatures);

centroids = msr.Center;
centroids = centroids(:,sortId); % to sort the centers based on Id
NRG = findNeighbors(diplbl, centroids, 35);

K = 2; %4
opts = statset('Display','final');
[idx,ctrs,sumdist] = kmeans(nRawFeatures,K,'dist','sqEuclidean', 'start', 'cluster',...
    'replicates',5,'EmptyAction','singleton','Options',opts);

% choose regino with low gray value which showthe germ mass.
% germCluster = 1;
% if ctrs(2,1) < ctrs(1,1)
%     germCluster = 2;
% end
% this is updated for the new versions
germCluster = 1;
if ctrs(2,1) > ctrs(1,1)
    germCluster = 2;
end
% overlayCluster(en_im,double(diplbl),idx , K, 'K-Means');
%%
dataCost = computeDataCost(ctrs, nRawFeatures);
dataCost = percentileStretch(dataCost);
% neighborhoodM = neighborGraph2Matrix(NRG);
[~, neighborhoodM] = computeDataCost(NRG, idx, 'RatioOfDifferentCluster');

[NumLbls NumSites] = size(dataCost);

% doing some experiment with graph cut
% NumSites = 4; % the number of regions
% NumLbls = 3; % the number of initial clusters
GH = GCO_Create(NumSites, NumLbls);
% the cost of assign label to each region
% dc = [0 9 2 0;      % Sites 1,4 prefer  label 1
%     3 0 3 3;      % Site  2   prefers label 2 (strongly)
%     5 9 0 5;];   % Site  3   prefers label 3
GCO_SetDataCost(GH, dataCost);
% define neighborhood
% neighbor = [0 1 0 0;     % Sites 1 and 2 connected with weight 1
%             1 0 1 0;     % Sites 2 and 3 connected with weight 1
%             0 1 0 2;     % Sites 3 and 4 connected with weight 2
%             0 0 2 0;];
GCO_SetNeighbors(GH, neighborhoodM); 
GCO_SetLabelCost(GH, ones(1,NumLbls));
energy = GCO_Expansion(GH,-1);
newLbl = GCO_GetLabeling(GH); % the labels

GCCluster = overlayCluster(en_im,double(diplbl),newLbl,K,{'save', outFile});

l2 = label(GCCluster==(3 - germCluster), 1, 100 , 0);
connectedGerm = connectGerms(GCCluster == germCluster, l2 > 0);
% overlayCluster(en_im, connectedGerm, [], 1, {'save', [outFile(1:end-4) 'GermMass01.tif']});
GMtuned = removeSmallGermRegion(lab(:,:,1),connectedGerm, [800, 7000, 0.6, 0.9]);
[lumen numL allLumen] = detectLumen(imgFile);
if numL > 0 
    Lumens = combineGMLumen(GMtuned,lumen, numL);
    LGM  = bwareaopen(or(GMtuned, Lumens), 20000);%or(GMtuned, Lumens);
    overlayCluster(en_im, LGM, [], 1, {'save', [outFile(1:end-4) 'GMLumen.tif']});
    overlayCluster(en_im, allLumen, [], 1, {'save', [outFile(1:end-4) 'GMLumenAll.tif']});
end
overlayCluster(en_im, GMtuned, [], 1, {'save', [outFile(1:end-4) 'TunedGermMass.tif']});

% l2 = label(GCCluster==2, 1, 0 , 100);
% dipshow(l2);
GCO_Delete(GH); % delete the object
% % % HL = GCO_ListHandles();
% % % for hl = 1: length(HL)
% % %     GCO_Delete(HL(hl));
% % % end

% detect boarder of tubules
edgLbl = label(and(edg, ~LGM), 2, 5); %~GMtuned
tmsr = measure(edgLbl, [], {'Feret', 'PodczeckShapes', 'MajorAxes'}, [],2);
m2o = msr2obj(edgLbl, tmsr, 'PodczeckShapes', 2);
% save([outFile(1:end-4) 'edg.mat'], 'edg', 'edgLbl','tmsr','m2o','LGM');
[approxLines lineMap flag] = tubuleBorder(edg, m2o, LGM);

sertoli = bwareaopen(im2mat(and(l2 > 0, ~LGM)), 20);
[closeGM sertoliLbl, sertoliCntr, ~, PAndN] = findClosestGM(sertoli, LGM, 100, lineMap, 0);
if flag 
showSerotli(0,closeGM, [outFile(1:end-4) 'connStared.tif'], 0, sertoliLbl, or(LGM,sertoli), en_im, PAndN);
else
    showSerotli(1,closeGM, [outFile(1:end-4) 'noedgeconnStared.tif'], 0, sertoliLbl, or(LGM,sertoli), en_im, PAndN);
end
%%
% now tune Sertoli cell cluster to find correct Sertoli
% 1- remove sertoli cells that compeletly fall into ger mass class
% % % sertoli = bwareaopen(im2mat(and(l2 > 0, ~connectedGerm)), 20);
% sertoliBoundary = imdilate(sertoli,strel('disk',1)) - sertoli;
% lblSertolliBoudary = bwlabel(sertoliBoundary, 8);
% % % [closeGM01 sertoliLbl] = findClosestGM(sertoli, connectedGerm, 100);
% % % showSerotli(closeGM01, [outFile(1:end-4) 'conn01.tif'], 0, sertoliLbl, or(connectedGerm,sertoli), en_im);

% sertoliBoundary = imdilate(sertoli,strel('disk',1)) - sertoli;
% lblSertolliBoudary = bwlabel(sertoliBoundary, 8);

% % % % sertoli = bwareaopen(im2mat(and(l2 > 0, ~LGM)), 20);
% % % % thre = getGrayCutOffbyPercentage(lab(:,:,1), .92);
% % % % [closeGM sertoliLbl, sertoliCntr, ~, PAndN] = findClosestGM(sertoli, LGM, 100,lineMap, 0);%, lab(:,:,1), thre);
% % % % showSerotli(closeGM, [outFile(1:end-4) 'connStared.tif'], 0, sertoliLbl, or(GMtuned,sertoli), en_im, PAndN);
% % ang = regionAngle(closeGM, sertoliCntr);
% % posAng = makeAnglePositive(ang,closeGM);
% % showSerotliAng(closeGM, posAng,sertoliCntr, or(GMtuned,sertoli) , en_im, PAndN)
%%%%%%%%%%%%save([imgFile(18:end-4) 'WS.mat']);

%%
% sertoli = bwareaopen(im2mat(and(l2 > 0, ~GMtuned)), 20);
% % sertoliBoundary = imdilate(sertoli,strel('disk',1)) - sertoli;
% % lblSertolliBoudary = bwlabel(sertoliBoundary, 8);
% thre = getGrayCutOffbyPercentage(lab(:,:,1), .92);
% [closeGM sertoliLbl, sertoliCntr, ~, PAndN] = findClosestGM(sertoli, GMtuned, 100);%, lab(:,:,1), thre);
% showSerotli(closeGM, [outFile(1:end-4) 'connStared.tif'], 0, sertoliLbl, or(GMtuned,sertoli), en_im, PAndN);


close all
end