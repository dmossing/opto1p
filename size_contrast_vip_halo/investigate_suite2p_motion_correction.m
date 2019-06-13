%%
result_fold = '/media/mossing/backup_0/data/suite2P/results/M10368/190410/2_3_4/suite2p/plane0';
raw_fold = '/media/mossing/backup_0/data/suite2P/raw/M10368/190410/2';
reg_fold = [result_fold '/reg_tif/'];
%%
cd(result_fold)
load('Fall.mat')
%%
cd(raw_fold)
fname = 'M10368_080_002_t00.tif';
info = imfinfo(fname);
num_images = 200; %numel(info);
img = imread(fname,1,'Info',info);
img = uint16(zeros([size(img) num_images]));
for k=1:num_images
    img(:,:,k) = imread(fname,1+4*(k-1),'Info',info);
end
%%
cd(reg_fold)
fname = 'file_chan000.tif';
info = imfinfo(fname);
num_images = numel(info);
img_pre_corr = imread(fname,1,'Info',info);
img_pre_corr = uint16(zeros([size(img_pre_corr) num_images]));
for k=1:num_images
    img_pre_corr(:,:,k) = imread(fname,k,'Info',info);
end
% %%
% imagesc(max(img,[],3))
%%
img_corr = img;
for k=1:num_images
    for j=1:size(ops.xblock,1)
        yinds = (1+ops.yblock(j,1):ops.yblock(j,2));
        xinds = (1+ops.xblock(j,1):ops.xblock(j,2));
        offset = round([ops.xoff1(k,j),ops.yoff1(k,j)]);
        img_corr(yinds,xinds,k) = circshift(img(yinds,xinds,k),offset);
    end
end
%%
for i=1:num_images, 
    subplot(1,3,1)
    imshow(img(:,:,i)),
    subplot(1,3,2)
    imshow(img_corr(:,:,i))
    subplot(1,3,3)
    imshow(img_pre_corr(:,:,i))
    pause
end
%%
for i=1:num_images
    red = img_corr(:,:,i);
    green = img_pre_corr(:,:,i);
    rg = cat(3,red,green,zeros(size(red)));
    imshow(rg)
    pause
end