%%
nplanes = 4;
fs = 7.75;
%%
roifile = cell(nplanes,1);
for i=1:nplanes
    %     roifile{i} = load(sprintf('M9835_250_004_ot_00%d.rois',i-1),'-mat');
    roifile{i} = load(sprintf('M9835_095_003_ot_00%d.rois',i-1),'-mat');
end
%%
for i=1:nplanes
    subplot(nplanes,1,i)
    hist(mean(roifile{i}.Neuropil),100)
end

%%
% hold on;
% [h1,x1] = hist(mn(idx==bigger),100);
% [h2,x2] = hist(mn(idx~=bigger),100);
% bar(x1,h1,'r');
% bar(x2,h2,'b');
% hold off;
%%
imagesc(roifile{1}.Data)
%%
% lkat = 102;
% t = (1:size(roifile{1}.Data,2))/fs;
% hold on;
% plot(t,roifile{1}.Data(lkat,:)')
% plot(t,artifact)
% plot(t,roifile{1}.Data(lkat,:)'-artifact(lkat,:)')
% hold off;
%%
% load('/home/mossing/modulation/visual_stim/190227/M9835/M9835_250_004.mat','result')
% load('/media/mossing/data_ssd/data/2P/190227/M9835/ot/M9835_250_004.mat','info
% load('/home/mossing/modulation/visual_stim/190301/M9835/M9835_175_002.mat','result')
% load('/media/mossing/data_ssd/data/2P/190301/M9835/ot/M9835_175_002.mat','info')
load('/home/mossing/modulation/visual_stim/190307/M9835/M9835_095_003.mat','result')
load('/media/mossing/backup_1/data/2P/190307/M9835/ot/M9835_095_003.mat','info')
%%
lights_on = result.gratingInfo.lightsOn(end,:);
%%
while find(diff(info.frame)<0,1)
    seam = find(diff(info.frame)<0,1);
    info.frame(seam+1:end) = info.frame(seam+1:end)+65536;
end
%%
offset = 1; % the first these many triggers are fake
for i=1:nplanes
    roiline = round(roifile{i}.ctr(1,:));
    if isfield(info,'rect')
        roiline = roiline + info.rect(1);
    end
    affected = zeros(size(roifile{i}.Data));
    for j=1:numel(lights_on)
        if lights_on(j)
            frames = 1+floor((-2+info.frame(offset+(j-1)*4+1))/nplanes):1+floor((-2+info.frame(offset+j*4))/nplanes);
            lines = [info.line(offset+(j-1)*4+1) info.line(offset+j*4)];
            lines(1) = lines(1) + mod(-2+info.frame(offset+(j-1)*4+1),nplanes)*512;
            lines(end) = lines(end) + mod(-2+info.frame(offset+j*4),nplanes)*512;
            affected(:,frames(2:end-1)) = 1;
            affected(:,frames(1)) = ((i-1)*512+roiline)>lines(1);
            affected(:,frames(end)) = ((i-1)*512+roiline)<lines(end);
        end
    end
    
    mn = mean(roifile{i}.Neuropil)';
    [idx,C] = kmeans(mn,2);
    [~,bigger] = max(C);
    artifact = abs(diff(C))*affected; %(idx==bigger)
    if ~isfield(roifile{i},'opto_stim_corrected')
        roifile{i}.Data = roifile{i}.Data-artifact; %repmat(artifact',size(roifile{i}.Data,1),1);
        roifile{i}.Neuropil = roifile{i}.Neuropil-artifact; %repmat(artifact',size(roifile{i}.Neuropil,1),1);
        roifile{i}.opto_stim_corrected = 1;
    end
end
%%
data = roifile{1}.Data; % - abs(diff(C))*affected;
neuropil = roifile{1}.Neuropil; % - abs(diff(C))*affected;
%
plot(data(1,:))
hold on;
plot(5000+1000*affected(1,:))
hold off
%%
% lkat = 101
% d = data(lkat,:)-abs(diff(C))*affected(lkat,:);
% n = neuropil(lkat,:)-abs(diff(C))*affected(lkat,:);
% plot(d-n)
%%
% for i=1:nplanes
%
% end
%%
for i=1:nplanes
    temp = roifile{i};
    temp.redratio = [roifile{i}.ROIdata.rois(:).redratio];
    save(sprintf('M9835_095_003_ot_00%d.rois',i-1),'-mat','-v7.3','-struct','temp');
end
%%
run_preprocessing_fold('.','/home/mossing/modulation/running/');
%%
% figure;
% ind = 8;
% af = affected(ind,:).*artifact(ind,:);
% % plot(affected(1,:).*artifact(1,:))
% % xlim([0 400])
% hold on;
% iplane = 3;
% plot(roifile{iplane}.Data(ind,:))
% plot(roifile{iplane}.Data(ind,:)-af)
% %%
% D = diff(roifile{1}.Data-affected.*artifact,[],2);
% figure;
% imagesc(D(:,1:500))
% %%
% figure;
% hist(D(:,37),100)
% hold on;
% hist(D(:,84),100,'facecolor','r')
% hold off
% %%
% figure;
% hist(D(:,36),100)
% hold on
% hist(D(:,37),100)
% hold off
% %%
% light_inds = find(lights_on);
% broke = [2 4 6 7 8 9 11 13]; % +1
% % broke = [1 2 3 5 6 7 9 10 12]; % +2
% j = light_inds(2); %broke(1));
% if lights_on(j)
%     frames = 1+floor((-2+info.frame(offset+(j-1)*4+1))/nplanes):1+floor((-2+info.frame(offset+j*4))/nplanes);
%     lines = [info.line(offset+(j-1)*4+1) info.line(offset+j*4)];
%     lines(1) = lines(1) + mod(-2+info.frame(offset+(j-1)*4+1),nplanes)*512;
%     lines(end) = lines(end) + mod(-2+info.frame(offset+j*4),nplanes)*512;
%     affected(:,frames(2:end-1)) = 1;
%     affected(:,frames(1)) = ((i-1)*512+roiline)>lines(1);
%     affected(:,frames(end)) = ((i-1)*512+roiline)<lines(end);
% end
% %%
% [min(((i-1)*512+roiline)>lines(1)) max(((i-1)*512+roiline)>lines(1))]
% [min(((i-1)*512+roiline)>lines(end)) max(((i-1)*512+roiline)>lines(end))]
% 
% %%
% frames