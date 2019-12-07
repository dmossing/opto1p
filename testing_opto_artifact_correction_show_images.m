%%
filename = 'file_chan050.tif';
N = 500;
im = zeros(508,654,N,'uint16');
for i=1:N
    im(:,:,i) = imread(filename,i);
end
%%
plot(squeeze(mean(mean(im))))
%%
for i=1:N
    imshow(im(:,:,i))
    pause(0.01)
end