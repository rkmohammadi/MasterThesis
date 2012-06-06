dipshow(lab(:,:,1));
out = dipprofile();


feildName = fieldnames(msr);
msrC1 = msr(class1Region);
msrC2 = msr(class2Region);

for fn =1 : length(feildName)
    t1 = getfield(msrC1, feildName{fn});
    t2 = getfield(msrC2, feildName{fn});
    switch size(t1,1)
        case 1
            h= scatterplot([t1; 1:length(t1)]');
            hold on
            scatterplot([t2; 1:length(t2)]',1,1,'+ red',h);
            title(feildName{fn});
            hold off
            
        case 2
            h= scatterplot(t1');
            hold on
            scatterplot(t2',1,1,'+ red',h);
            title(feildName{fn});
            hold off
        otherwise
            fprintf(1,'the size of %s is %d,%d \n',feildName{fn}, size(t1));
    end
end

%% show the selected region
for i =1 : numel(images)
    imgFile = images{i};
    load([imgFile(1:end-3) 'mat']);
    
    figure(91),
    imshow(en_im);
    hold on; 
    mask = zeros(size(lbl_bin_FRST));
    for j =1 : length(class1Region)
        if class1Region(j) ~= 0
            mask = mask + (lbl_bin_FRST== class1Region(j));
        end
    end
    mask = mask* 255;
    colorMask = uint8(cat(3, mask, mask, mask));
    hold on
    handel = imshow(colorMask); title(['regions on top of the image', imgFile]);
    set(handel,'AlphaData',0.4);
    hold off
end