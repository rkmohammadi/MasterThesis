function [lumen numL lsBigLumen]= detectLumen(I)
% function [lumen numL lsBigLumen]= detectLumen(I)
% This function get a grey level image and select bright regions.
% the method is described in 4.2 Detect Lumen.
% Input:
%   I, a grey level image. For this project I have used the L-channel of
%   Lab color space.
% Output:
%   lumen, all bright region that fulfill the size and convexity
%   constraint.
%   numL, the number of accepted bright regions that are in lumen.
%   lsBigLumen, all birght region after level set method with checking the
%   size and convexity of the regions.
% 
% 
% Author: A.Rahim Kadkhodamohammadi (r.k.mohammadi@gmail.com)
% 01 June 2012 CBA, Uppsala University
%--------------------------------------------------------------

initThre = 0.8;

if ischar(I)
    lab = RGB2Lab(imread(I));
    l_channel = lab(:,:,1);
else
    l_channel = I;
end
sz = size(l_channel);

w = 5;
sigma = [3, .2];
l_channel = bfilter2(l_channel./max(l_channel(:)),w,sigma);

level = getGrayCutOffbyPercentage(l_channel, .6);
l01 = (l_channel > level);

% % %  since we use levelset we are interested to have smaller regions
% % % openL01 = opening(l01,15);
% % % bigLumen = label(openL01,2, 8000,0);
% % % bigLumen = closing(bigLumen, 10); % bigLumen>0
%lblLumen = bigLumen; % label(bigLumen);

% rpalcement code for levelset
erodedL01 = erosion(l01,10);
openLumen = opening(erodedL01);
bigLumen = label(openLumen,2,8000,0);
[~, lsBigLumen] = lumen01(l_channel, bigLumen>0);
lsBigLumen = label(logical(lsBigLumen));
msr = measure(lsBigLumen, [], {'Size','Convexity'},[],2);
ind = [];
% the ratio of changing threshold value for choosing lumen will be
% changed by alpha value. the initial value of threshold is in
% initThre. it will be tuned based on the size of lumen and alpha.

alpha = 0.041; % .33
for jj = 1: size(msr,1)
    thre = initThre;
    if msr(jj).Size< 150000
        thre = thre - msr(jj).Size*alpha/25000;
    else
        thre = 0.6;
    end
    if msr(jj).Convexity >= thre
        ind = [ind, jj];
    end
end

lumen= zeros(sz(1:2));
numL = numel(ind);
for r =1 : numL
    lumen(lsBigLumen == ind(r)) = r;
end
end