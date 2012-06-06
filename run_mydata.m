%% create ucm
addlibs();

folder = '/home/rahim/work/Dropbox/images/GATA-4/';
files = dir([ folder '*.TIF']);
for i = numel(files) : -1 : 1   
    imgFile = [folder files(i).name];
    outFile = [folder 'gPb_' files(i).name(1:end-3) 'mat'];
    % preprocessing
    Im = imread(imgFile);
    doubleIm = double(Im) ./255;
    w = 5;
    sigma = [3, .2];
    en_im = bfilter2(doubleIm, w, sigma);

    gPb_orient = globalPb(en_im, outFile, 0.3);
    
    ucm = contours2ucm(gPb_orient, 'imageSize');
    imwrite(ucm,[folder 'ucm_' files(i).name(1:end-3) 'bmp'], 'bmp');
end

%%

for i = numel(files)-1 : -5 : 1   
    imgFile = [folder files(i).name];
    Im = imread(imgFile);
    sz = size(Im);
    % preprocessing
    doubleIm = double(Im) ./255;
    w = 5;
    sigma = [3, .2];
    en_im = bfilter2(doubleIm, w, sigma);
    
    % fast radial symmetry of s plane in HSV
    hsv_im = rgb2hsv(en_im);
    
    % convert s channel to uint8
    s_channel = percentileStretch(hsv_im(:,:,2), [0, 1]);
    
    gaussK = fspecial('gaussian',[7 7],2);
    s_channel = imfilter(s_channel, gaussK,'same');
    
    [S, So] = fastradial(s_channel, [3 4 5], 1, .2);
    Bin_s = S > 1; 
    dipshow(Bin_s);
    figure, imshow(Im);
    InterestPoint = uint8(cat(3, (Bin_s)*255 ,(Bin_s)*255,(Bin_s)*255));
    hold on
    distIm = imshow(InterestPoint); title(['regions on top of the image', imgFile]);
    set(distIm,'AlphaData',0.4);
    
%     pause
end
%%
for i = numel(files) : -1 : 1   
    imgFile = [folder files(i).name];
    % preprocessing
    Im = imread(imgFile);
    lab = RGB2Lab(Im);
    w = 5;
    sigma = [3, .2];
    rowl = bfilter2(lab(:,:,1),w,sigma);
    level = getGrayCutOffbyPercentage(rowl, .6);
    lumen = (rowl > level);
    openLumen = opening(lumen);
    bigLumen = label(openLumen,2, 5000,0);
    bigLumen = closing(bigLumen>0, 10);
    lblLumen = label(bigLumen);
    msr = measure(lblLumen, [], 'Convexity',[],2);
    convex = msr.Convexity;
    mId = msr.ID;
    [~, sortm] = sort(mId);
    convex = convex(sortm);
    ind = find(convex>0.7);
    sz = size(Im);
    clusterdRgn= zeros(sz(1:2));
    outFile = ['../../tmpResult/' files(fi).name(1:end-3) 'lumen01.tif'];
    overlayCluster(Im, lblLumen >0 , [], 1, {'save', outFile});
    for r =1 : numel(ind)
                 clusterdRgn(lblLumen == ind(r)) = 1;
    end
    outFile = ['../../tmpResult/' files(fi).name(1:end-3) 'lumen02.tif'];
    overlayCluster(Im, clusterdRgn, [], 1, {'save', outFile});
    close all
end
