%%
cd('/media/mossing/data_ssd/data/2P/190320/M10365')
im1 = sbxread('M10365_125_001',0,10);
im2 = sbxread('M10365_125_002',0,10);
%%
figure
for i=1:10
subplot(1,2,1)
imshow(squeeze(im1(1,:,:,i)))
subplot(1,2,2)
imshow(squeeze(im2(1,:,:,i)))
pause
end