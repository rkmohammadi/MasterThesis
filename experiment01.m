%%
%% smoothing with gaussian kernel
gaussK = fspecial('gaussian',[3 3],1);
Ig = imfilter(a, gaussK,'same');
figure, imshow(Ig), title('smoothed saturation with Gaussian');

 figure, imshow(im2bw(Ig, graythresh(Ig))),