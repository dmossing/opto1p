% [offset_image,flowed_image] = function align_sub_images(im1,im2,ops)
% ops should have fields xblock, yblock, as in the convention used by
% suite2p
% computes a flow field offset_image, s.t.
% flowed_image = apply_flowfield(offset_image,im2) maximally resembles im1. 
im1 = ops.meanImg;
im2 = green_mean;
offset = gui_optimize_image_offset(im1,im2);
%%
figure
rg = combine_rg(im1,circshift(im2,[v,u]),1);
imshow(rg)
%%
redundancy = zeros(size(green_mean));
for i=1:size(ops.yblock,1)
    yrg = ops.yblock(i,1)+1:ops.yblock(i,2);
    xrg = ops.xblock(i,1)+1:ops.xblock(i,2);
    redundancy(yrg,xrg) = redundancy(yrg,xrg)+1;
end
%%
offsets = zeros(size(ops.yblock,1),2);
A = ops.meanImg;
B = circshift(green_mean,offset);
offset_image = zeros([size(green_mean) 2]);
% registered_image = zeros(size(green_mean));
for i=1:size(ops.yblock,1)
    yrg = ops.yblock(i,1)+1:ops.yblock(i,2);
    xrg = ops.xblock(i,1)+1:ops.xblock(i,2);
    [u,v] = fftalign_dpm(A(yrg,xrg),B(yrg,xrg));
    offsets(i,:) = [u v];
    for idim=1:2
        offset_image(yrg,xrg,idim) = offset_image(yrg,xrg,idim) - offset(idim) + offsets(i,idim);
    end
%     registered_image(yrg,xrg) = registered_image(yrg,xrg) + circshift(B(yrg,xrg),-offsets(i,:));
%     offsets(i,:) = gui_optimize_image_offset(ops.meanImg(yrg,xrg),green_mean(yrg,xrg),init_offset);
%     init_offset = offsets(i,:);
end
offset_image = offset_image./repmat(redundancy,[1 1 2]);
% registered_image = registered_image./redundancy;
% %%
% offsetbox = permute(reshape(offsets,[8 6 2]),[2 1 3]);
% %%
% figure
% subplot(1,2,1)
% imagesc(offset_image(:,:,1))
% subplot(1,2,2)
% imagesc(offset_image(:,:,2))
%%
flowed_image = apply_flowfield(green_mean,offset_image);
% flowed_red = apply_flowfield(red_mean,offset_image);
% %%
% figure
% % subplot(1,3,1)
% % rg = combine_rg(ops.meanImg,registered_image,1);
% % imshow(rg)
% 
% % subplot(1,3,2)
% rg = combine_rg(flowed_red,ops.meanImg,1);
% imshow(rg)
% 
% % subplot(1,3,3)
% % rg = combine_rg(flowed_image,registered_image,1);
% % imshow(rg)