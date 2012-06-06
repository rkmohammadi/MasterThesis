folder = '/home/rahim/work/Dropbox/images/GATA-4/';
files = dir([ folder '*.TIF']);
for i = numel(files) : -1 : 1   
    imgFile = [folder files(i).name];
    % preprocessing
    Im = imread(imgFile);
    lab = RGB2Lab(Im);
    rowl = lab(:,:,1);
    w = 5;
    sigma = [3, .2];
    rowl = bfilter2(rowl./max(rowl(:)),w,sigma);
    level = getGrayCutOffbyPercentage(rowl, .6);
    lumen = (rowl > level);
    openLumen = opening(lumen,15);
    bigLumen = label(openLumen,2, 8000,0);
    bigLumen = closing(bigLumen, 10); % bigLumen>0 
    lblLumen = bigLumen; % label(bigLumen);
    msr = measure(lblLumen, [], {'Size','Convexity'},[],2);
    convex = msr.Convexity;
    %mId = msr.ID;
    %[~, sortm] = sort(mId);
    %convex = convex(sortm);
    ind = [];
    alpha = 0.033;
    for jj = 1: size(msr,1)
        thre = 0.8;
        if msr(jj).Size< 150000
            thre = thre - msr(jj).Size*alpha/25000;
        else
            thre = 0.6;
        end
        if msr(jj).Convexity > thre
            ind = [ind, jj];
        end
    end
    % = find(convex>0.7);
    sz = size(Im);
    clusterdRgn= zeros(sz(1:2));
    outFile = ['../../tmpResult/' files(i).name(1:end-3) 'lumen01.tif'];
    overlayCluster(Im, lblLumen >0 , [], 1, {'save', outFile});
    for r =1 : numel(ind)
                 clusterdRgn(lblLumen == ind(r)) = 1;
    end
    outFile = ['../../tmpResult/' files(i).name(1:end-3) 'lumen02.tif'];
    overlayCluster(Im, clusterdRgn, [], 1, {'save', outFile});
    close all
end

%%

%Finding edges of a color image
%Authors : Jeny Rajan, Chandrashekar P.S 
%Usage edgecolor('abc.jpg');

%%function R=edgecolor(nm);
% % img=imread(nm);
[x y z]=size(en_im);
if z==1
    rslt=edge(img,'canny');
elseif z==3
    img1=rgb2ycbcr(en_im);
    dx1=edge(img1(:,:,1),'canny');
    dx1=(dx1*255);
    img2(:,:,1)=dx1;
    img2(:,:,2)=img1(:,:,2);
    img2(:,:,3)=img1(:,:,3);
    rslt=ycbcr2rgb(uint8(img2));
end
R=rslt;

%%
    enLine = dip_pathopening(lab(:,:,1), ~GMtuned, 50,0,0);
    BW = and(edge(im2mat(enLine), 'canny', [], 1 ), ~GMtuned);
    enLine1 = dip_pathopening(II(:,:,1), [], 50,0,0);
    enLine2 = dip_pathopening(II(:,:,2), [], 50,0,0);
    enLine3 = dip_pathopening(II(:,:,3), [], 50,0,0);
    BW1 = and(edge(im2mat(enLine1), 'canny', [], 1 ), ~GMtuned);
    BW2 = and(edge(im2mat(enLine2), 'canny', [], 1 ), ~GMtuned);
    BW3 = and(edge(im2mat(enLine3), 'canny', [], 1 ), ~GMtuned);
    
    
    dipshow(cat(3, enLine, enLine1,enLine2,enLine3), 'all')
    
     dipshow(cat(3, BW, BW1,BW2,BW3), 'all')
%%

% % I  = imread('circuit.tif');
% % rotI = imrotate(I,33,'crop');
% % BW = edge(rotI,'canny');
BW = and(edge(im2mat(enLine), 'canny', [], 1 ), ~GMtuned);
[H,T,R] = hough(BW,'RhoResolution',.5,'Theta',-90:.5:89.5);
imshow(H,[],'XData',T,'YData',R,...
            'InitialMagnification','fit');
xlabel('\theta'), ylabel('\rho');
axis on, axis normal, hold on;
P  = houghpeaks(H,100,'threshold',ceil(0.4*max(H(:))));%, 'NhoodSize',[15,7]);
x = T(P(:,2)); y = R(P(:,1));
plot(x,y,'s','color','white');
% Find lines and plot them
lines = houghlinesadopted(BW,T,R,P,'FillGap',30,'MinLength',110,'tol',[6 0]);
h = figure, imshow(en_im), hold on
max_len = 0;
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

   % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > max_len)
      max_len = len;
      xy_long = xy;
   end
end
%%
for ii = 2 : numel(files)
    load(['../../tmpResult/' files(ii).name(1:end-4) 'edg.mat']);
    dipshow(1,edg, 'all');
    dipshow(2,and(m2o> 0, m2o<0.08))
    pause; 
end

%%
% function [ O ] = hough_transform( fPath )
%   fPath: the path to the image
%   this program implements the Standard Hough Transform 
%   to detect lines
    rhoStep=1;
    thetaStep=1;
    rhoDiffThresholdForLine=rhoStep/8;
    I=imread(fPath, 'gif');
    %find the edge of the image
    BW=edge(I,'canny');
    %imshow(BW);
    %hough transform
    %define the accumulator range
    rho=1:rhoStep:sqrt((size(BW,1))^2 + (size(BW,2))^2);
    theta=0:thetaStep:180-thetaStep;
    accu=zeros(length(rho), length(theta));
    %get the pixel indices that contains a point
    [rowInd, colInd]=find(BW);
    %for each point, plot all the lines (sampled) pass through it
    %at theta-rho plane
    for li=1:1:length(rowInd)
        for lk=1:1:length(theta)
            ltheta=theta(lk)*pi/180;
            lrho=colInd(li)*cos(ltheta) + rowInd(li)*sin(ltheta);
            %binning the lrho value
            diffs=abs(lrho-rho);
            %we only increase the count of most similar ones
            %introducing a threshold instead choosing the
            %min
            minDiff=min(diffs);
            if (minDiff<rhoDiffThresholdForLine)
               minDiffInd=find(diffs==minDiff);
               for lm=1:1:length(minDiffInd)
                   accu(minDiffInd(lm),lk) = accu(minDiffInd(lm),lk) + 1;
               end
            end
        end
    end
    %find local maxima 
    accuBMax=imregionalmax(accu);
    [rho_candi, theta_candi]=find(accuBMax==1);
    %find the points in theta-rho plane that has count more than
    %threshold
    linePoints=0;
    %get a list of lines detected with their rho and theta values
    rhoLines=[];
    thetaLines=[];
    for li=1:1:length(rho_candi)
        l_accu=accu(rho_candi(li), theta_candi(li));
        if (l_accu<=0)
            %do nothing
        elseif (l_accu > 25)
            linePoints=linePoints+1;
            rhoLines=[rhoLines;rho(rho_candi(li))];
            thetaLines=[thetaLines;theta(theta_candi(li))];
        end
    end
    end

    
