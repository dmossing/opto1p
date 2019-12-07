%%
filename = 'file_chan050.tif';
N = 500;
im = zeros(508,654,N,'uint16');
for i=1:N
    im(:,:,i) = imread(filename,i);
end
%%
N = 500;
cd('/home/mossing/recurrence/2P/191017/M0276')
im = sbxreadpacked('M0276_120_003',1000,N);
%%
im = im(:,:,1:4:end);
%%
plot(squeeze(mean(mean(im))))
%%
for i=1:N/4
    imshow(im(:,:,i))
    pause(0.01)
end
%%
imshow(im(:,:,40))
%%
ii = reshape(im(:,151:end-150,:),[],size(im,3));
%%
scatter(ii(:,38),ii(:,40),1)
%%
t1 = 58;
t2 = 59;
figure
subplot(1,2,1)
imshow(im(:,151:end-150,t1))
subplot(1,2,2)
imshow(im(:,151:end-150,t2))
%%
imshow(cat(3,im(:,:,t1),im(:,:,t2),0*im(:,:,1)))
%%
figure
hold on
scatter(ii(:,t1),double(ii(:,t2))-mean(ii(:,t2))+mean(ii(:,t1)),1)
% scatter(ii(:,t1),double(ii(:,t2)),1)
hold off
% plot([0 3e4],[0 -3e4])
% hold off