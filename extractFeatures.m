function [cells, features]= extractFeatures(I, L, stats, Image) %center, boundingBox, pixelList)
% it is the function I got from Azadeh but I have not used

% I imput image in rgb
% L labeled image 
% 

% I1 =double(imread('M1x10ex1.tif'));
% h = fspecial('gaussian', 7, 1);
% I2=imfilter(I1,h,'same');
% I=imresize(I2,0.5);
% 
sz=size(I);

% r=I(:,:,1);
% a(:,1)=r(:);
% g=I(:,:,2);
% a(:,2)=g(:);
% b=I(:,:,3);
% a(:,3)=b(:); 

% d=joinchannels('rgb',I);



K = 2;
%% Kmenas 2 clsuetr 
% X = a; K = 2;
% opts = statset('Display','final');
% [idx,ctrs,sumdist] = kmeans(X,K,'dist','sqEuclidean',...
% 'replicates',10,'EmptyAction','singleton','Options',opts);
% out=reshape(idx,sz(1),sz(2));
% figure,imagesc(out)
% m=mean(ctrs,2)
% o=find(m==min(m));
% tt=zeros(sz(1),sz(2));
% tt(out==o)=1;
% L=label(tt>0,1,20,500)

% L = label(BWI, 1, 0, 0);


% finding features
%  msr1 = measure(L,d{1},({ 'StdDev', 'perimeter','mean','size','center'}));
%  msr2 = measure(L,d{2},({ 'StdDev', 'perimeter','mean'}));
%  msr3 = measure(L,d{3},({ 'StdDev', 'perimeter','mean'}));
%  msr=[msr1.mean;msr2.mean;msr3.mean;msr1.stdDev;msr2.stdDev;msr3.stdDev];
%  msr_mean=mean(msr);
 
 
%  center = round(msr1.center);
 num = numel(stats); % the number of regions
 cells = cell(num,1);
 obj = cell(num, 1);
 ent = zeros(num, sz(3));
 msr = zeros(num, 2*sz(3));
%% 
 for i = 1 : num
     boundingBox = round(stats(i).BoundingBox);
     lblBoundingB = [];
     lblBoundingB = L(boundingBox(2): min(boundingBox(2)+boundingBox(4), sz(1)), ...
         boundingBox(1): min(boundingBox(1)+boundingBox(3),sz(2)),:);
     lblBoundingB = (lblBoundingB==i);
     lblBoundingB = repmat(lblBoundingB,[1,1,sz(3)]);
     o1 = []; 
     o1 = I(boundingBox(2): min(boundingBox(2)+boundingBox(4), sz(1)), ...
         boundingBox(1): min(boundingBox(1)+boundingBox(3), sz(2)),:);
     cells{i} = uint8(double(o1).*(lblBoundingB));
     
     PL = stats(i).PixelList;
     obj1 = zeros(length(PL), sz(3));
     for k = 1 : sz(3)
         for j=1 : length(PL)
             obj1(j ,k) =I(PL(j,2), PL(j,1), k);
         end
         msr(i,k) = mean(obj1(:,k));
         msr(i,sz(3)+k) = std(obj1(:,k));
         
         % compute intensity entropy for each channel
         vr=round(obj1(:,k));
         [pdf_obj1,bins]=hist(vr,[min(vr):5:max(vr)]);
         
         pdf_obj=pdf_obj1/sum(pdf_obj1);
         in=find(pdf_obj==0);
         if ~isempty(in)
             pdf_obj(in)=[];
         end
         ent(i,k)=-sum(pdf_obj.*log2(pdf_obj));
     end
     obj{i} = obj1;
 end
 
     
%  for i=1:max(max(L)) 
%      [x,y]=find(double(L)==i);
%      obj1=[];
%      o1=zeros(sz(1),sz(2),3);
%      for j=1 :length(x)         
%       obj1(j,1)=I(x(j),y(j),1);
%       obj1(j,2)=I(x(j),y(j),2);
%       obj1(j,3)=I(x(j),y(j),3);
%       o1(x(j),y(j),:)=I(x(j),y(j),:);
%      end
%          
%      obj{i}=obj1;
%      indx1=max(center(2,i)-30,1);
%      indx2=min(center(2,i)+30,sz(1));
%      indy1=max(center(1,i)-30,1);
%      indy2=min(center(1,i)+30,sz(2));     
%      cells{i}=o1([indx1:indx2],[indy1:indy2],:);  
%  end

 
 % compute entropy of each region
%  for i=1:max(max(L))
%      for j=1:3         
%            v=obj{i}(:,j);
%            vr=round(v);
%            [pdf_obj1,bins]=hist(vr,[min(vr):5:max(vr)]);
%            
%            pdf_obj=pdf_obj1/sum(pdf_obj1);
%            in=find(pdf_obj==0);
%            if ~isempty(in) %length(in)~=0
%              pdf_obj(in)=[];
%            end 
%            
%            ent(i,j)=-sum(pdf_obj.*log2(pdf_obj));
%            
%      end
%  end
 
 
%% 
% cooccurance
numG=3;
of1=[zeros(numG,1),(1:numG)'];
of2=[-(1:numG)',(1:numG)'];
of3=[-(1:numG)', zeros(numG,1)];
of4=[-(1:numG)',-(1:numG)'];
offset=[of1;of2;of3;of4];
szof=size(offset);
numlevel=8;
energy = zeros(length(obj),sz(3));
hmgnty = zeros(length(obj),sz(3));
cor = zeros(length(obj),sz(3));
entropy= zeros(length(obj),sz(3));

warningId = 'images:graycomatrix:scaledImageContainsNan';
warning('off', warningId);
for i=1:length(obj) 

    for j=1:sz(3)

    input=cells{i}(:,:,j);
    input(input==0)=NaN;
    glcm = graycomatrix(input,'NumLevels',numlevel,'Offset',offset,'GrayLimits',[]);
    stats=graycoprops(glcm,'Energy','Homogeneity','Correlation');    
    energy(i,j)=mean(stats.Energy);
    hmgnty(i,j)=mean(stats.Homogeneity);
    TF = isnan(stats.Correlation);
    stats.Correlation(TF)=[];
    cor(i,j)=mean(stats.Correlation);

    for jj=1:szof(1)
        jointProb=glcm(:,:,jj);
        if (sum(sum(jointProb))~=0)
          JP=jointProb/sum(sum(jointProb));
        else
            JP=0;
        end
        
       
        [indx,indy]=size(JP);
        x=repmat((1:indx)',indy,1);
        y=repmat((1:indy),indx,1);
        y=y(:);
        JP=JP(:);          
        
        DM1=sum((x-y).*2.*JP);
        in=find(JP==0);
           if ~isempty(in) %length(in)~=0
             JP(in)=[];
           end 
           entropy1=-sum(JP.*log2(JP));
      
    end
   
    DM(i,j)=mean(DM1);
    entropy(i,j)=mean(entropy1);  
   
  
    end
end

 
%% 
features=[msr,ent,energy,hmgnty,entropy,DM,cor];
[n,p]=size(features);
features = (features - ones(n,1)*mean(features)) ./ (ones(n,1)*std(features));

[pc score latent]=princomp(features);
rvar1=cumsum(latent)/sum(latent);
figure,plot(rvar1);
grid on
legend('retaind variance');


X=score(:,1:10);
[n,p]=size(X);
opts = statset('Display','final');
k=[1,2,3,4,5,6,7,8,9,10];
for j=1:10
    for i=1:10
         [idx,ctrs,sumdist] = kmeans(X,k(j),'dist','sqEuclidean',...
         'replicates',1,'EmptyAction','singleton','Options',opts);
         SW(i)=mean(sumdist);
    end
   
SWk(j,:)=(SW);
% std_SW(j)=std(SW);
end
SWkn=SWk./repmat(SWk(1,:),10,1);
meanSW=mean(SWkn');
stdSW=std(SWkn');
figure,errorbar(k,meanSW,stdSW);
grid on



%%
 K = 4; %4
opts = statset('Display','final');
[idx,ctrs,sumdist] = kmeans(X,K,'dist','sqEuclidean', 'start', 'cluster',...
'replicates',5,'EmptyAction','singleton','Options',opts);

% for j=1:K
% l=find(idx==j);
% y=zeros(sz(1),sz(2));
% for i=1:length(L) %tt
%     y(L==l(i))=1;
% end
% figure,imshow(y)
% end
% figure,imshow(uint8(I));


sgm = repmat(eye(p),[1 1 K]);
priorp=ones(1,K)/K;
S=struct('mu',ctrs,'Sigma',sgm,'PComponents',priorp);
    
opts = statset('Display','final','MaxIter',10000,'TolFun',1e-20);
object = gmdistribution.fit(X,K,'Options',opts,'Replicates',1,'Regularize',1e-5,'Start',S);
idx = cluster(object,X);


for j=1:K   
l=find(idx==j);
y=zeros(sz(1),sz(2));
for i=1:length(l)
    y(L==l(i))=1;
end
figure,imshow(Image)
RGBDist = uint8(cat(3, zeros(sz(1:2)),(y)*255 ,zeros(sz(1:2))));
hold on
distIm = imshow(RGBDist); title(['Cluster ', int2str(j)]);
set(distIm,'AlphaData',0.3);
end
% figure,imshow(uint8(I));
% 
% hold on
% contour(y)
% hold off 
%%


ctr=object.mu;
mc=mean(ctr,2);
[ss,ss1]=sort(mc);
%ss1(2) is sertoli center
for i=1:K
    dist(i,:)=sqrt(sum((X-repmat(ctr(i,:),n,1)).^2,2));
end
expo=2;
tmp = dist.^(-2/(expo-1)); 
mf = tmp./(ones(K, 1)*sum(tmp));  
fuzzy=zeros(sz(1),sz(2));    

for i=1:max(max(L))
    fuzzy(L==i)=mf(ss1(2),i);    
end

figure,imshow(fuzzy,[]);

%membeship function higher than 0.5
tt=find(mf(ss1(2),:)>0.5);
y1=zeros(sz(1),sz(2));
for i=1:length(tt)
    y1(L==tt(i))=1;
end
figure,imshow(y1);
figure,imshow(uint8(I));
hold on
contour(y1);
