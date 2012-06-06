function [phi,white]=lumen01(I, Lerod)
% function [phi,white]=lumen01(I, Lerod)
% a warper function to set up levle set,
%  Input:
%   I, a gray levle image.
%   Lerod, binary image to initialize level-set.
%  Output:
%   phi, the output of level-set for growing the regions
%   white, the final output of level-set for the initial regions.
%
% Author: Azadeh
%-------------------------------------------------------------------------

% P=path;
% I=readim('R010007.TIF')



im = double(I);
sz = size(im);
%% added
% rgb=joinchannels('rgb',I);
% lab_I=colorspace(rgb,'lab');


%%
% % % % sz = size(im)
% % % % data = reshape(im, [prod(sz(1:2)) 2]);
% % % % [idx,c]=kmean3(data);
% % % % LMean = reshape(idx,sz(1),sz(2)); 
% % % % [maxind,temp]=find(c(:,1)==max(c(:,1))); 
% % % % %imagesc(LMean)
% % % % 
% % % % 
% % % % Ldip=dip_image(LMean);
% % % % lum=Ldip==maxind;
% % % % luml=label(lum,1,100);
% % % % L1=imfill(double(luml));
% % % % L=ones(sz(1),sz(2));
% % % % L(L1==0)=0;
% % % % se=strel('disk',1);  %%%%%one parameter here how much enrode
% % % % Lerod=imerode(L,se);



%figure,imshow(I)
% hold on;
% contour(L, 'k','LineWidth',1);
% figure,imshow(I)
% hold off
% path(P,'/home/azadeh/Documents/Personal/My-project/levelset/DRLSE_v0');
c0=2;
 initialLSF = c0*ones(sz(1),sz(2));
 initialLSF(Lerod==1)=-c0;

%  sigma=1;
%  Is=double(gaussf(I(:,:,2),sigma));
Is=I;%I(:,:,1); % 2
%% added
%Is=double(lab_I{2});
%%
Img=imresize(Is,0.5);
initialLSF=imresize(initialLSF,0.5);


%%%%edge indicator
% sigma=.8;    % scale parameter in Gaussian kernel
% G=fspecial('gaussian',15,sigma); % Caussian kernel
% Img_smooth=conv2(Img,G,'same');  % smooth image by Gaussiin convolution
% [Ix,Iy]=gradient(Img);

Ix=dx(Is);
Iy=dy(Is);
f=double(Ix).^2+double(Iy).^2;
g=1./(1+f);  % edge indicator function.
g=imresize(g,0.5);
% figure,imshow(f,[])
% dip_image(f)
% 
% Ix1=double(dx(I(:,:,1)));
% Ix2=double(dx(I(:,:,2)));
% Ix3=double(dx(I(:,:,3)));
% 
% Iy1=double(dx(I(:,:,1)));
% Iy2=double(dx(I(:,:,2)));
% Iy3=double(dx(I(:,:,3)));
% 
% g11=Ix1.^2+Ix2.^2+Ix3.^2;
% g12=Ix1.*Iy1+Ix2.*Iy2+Ix3.*Iy3;
% g22=Iy1.^2+Iy2.^2+Iy3.^2;
% 
% delta=(g11-g22).^2+4*(g12.^2);
% lamda_p=(g11+g22+sqrt(delta))/2;
% lamda_n=(g11+g22-sqrt(delta))/2;
% f=sqrt(lamda_p-lamda_n);
% g=1./(1+f);  % edge indicator function.
% %g=imresize(g,0.5);
%  figure,imshow(f,[]) 
%  dip_image(f)

%tic
itr=20; %%%one parameter here
[phi]=demo_1(Img,initialLSF,5,g); %6 parameter here.
%toc

%  figure,imshow(I)
 phi_re=imresize(phi,2);
%  hold on;  contour(phi_re, [0,0], 'k','LineWidth',1);
 phi_dip=dip_image(phi_re);
 white=double(phi_dip<0);









