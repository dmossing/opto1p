% %%
%%
endat = 10000;
offset0 = 11;
% tv = zeros(21,1);
% for offset=1:21
%     tv(offset) = sum(sum(abs(diff(roifile{i}.Data(:,end-endat-offset0:end-offset0)-artifact(:,end-endat-offset:end-offset),[],2))));
% end
% figure;
% plot(tv)
%%


%%
[~,minind] = min(tv1);
loffset1 = minind-noffset-1;
[~,minind] = min(tv2);
loffset2 = minind-noffset-1;

%%
figure;
% plot(-noffset:noffset,tv)
plot(tv1)
hold on;
plot(tv2)
hold off;

% %

artifact_size = max(abs(artifact_cell{i}(:)))*(2*(max(artifact_cell{i}(:))>0)-1);
affected = artifact_cell{i}/artifact_size;

[ii,jj] = find(diff(affected,[],2)>0);
T = 31;
Toffset = 10;

figure;
temp_data = roifile{i}.Data - 1*artifact_cell{i};
temp_neuropil = roifile{i}.Neuropil -1*artifact_cell{i};
for iii=1:30
    data_avg = zeros(1,T);
    neuropil_avg = zeros(1,T);
    af_avg = zeros(1,T);
    for j=1:numel(ii)
        if ii(j)==iii
            data_avg = data_avg + temp_data(ii(j),jj(j)-Toffset:jj(j)+T-Toffset-1)/numel(ii);
            neuropil_avg = neuropil_avg + temp_neuropil(ii(j),jj(j)-Toffset:jj(j)+T-Toffset-1)/numel(ii);
            af_avg = af_avg + affected(ii(j),jj(j)-Toffset:jj(j)+T-Toffset-1)/numel(ii);
        end
    end
    subplot(3,10,iii)
    plot(data_avg)
    hold on;
    plot(neuropil_avg)
    plot(af_avg*1000);
    hold off
end
%%

[ii,jj] = find(diff(control,[],2)>0);
T = 31;
Toffset = 10;
data_avg = zeros(1,T);
neuropil_avg = zeros(1,T);
af_avg = zeros(1,T);
figure;
for iii=1:10
    for j=1:numel(ii)
        if ii(j)==iii
            data_avg = data_avg + roifile{i}.Data(ii(j),jj(j)-Toffset:jj(j)+T-Toffset-1)/numel(ii);
            neuropil_avg = neuropil_avg + roifile{i}.Neuropil(ii(j),jj(j)-Toffset:jj(j)+T-Toffset-1)/numel(ii);
            af_avg = af_avg + affected(ii(j),jj(j)-Toffset:jj(j)+T-Toffset-1)/numel(ii);
        end
    end
    subplot(1,10,iii)
    plot(data_avg)
    hold on;
    plot(neuropil_avg)
    plot(af_avg*1000);
    hold off
end

%%
iroi = 5;
figure
hold on
% plot(roifile{1}.Data(iroi,:)); plot(artifact_cell{i}(iroi,:))
plot(temp_data(iroi,:)); plot(artifact_cell{i}(iroi,:))
hold off

%%
[mskcdf,iplane] = compute_mskcdf(roifile,info);
[data,~] = append_roifile_parts(roifile,'Data',false);
[neuropil,~] = append_roifile_parts(roifile,'Neuropil',false);
[ctr,iplane] = append_roifile_parts(roifile,'ctr',true);
roiline = round(ctr(1,:));
roiline = roiline(:);
if isfield(info,'rect')
    roiline = roiline + info.rect(1);
end

artifact = compute_artifact_mskcdf(mskcdf,iplane,neuropil,info,lights_on,loffset1,loffset2);
% artifact = compute_artifact_mskcdf(mskcdf,iplane,neuropil,info,lights_on,loffset1,loffset2+0);
% artifact = compute_artifact_fn(roiline,iplane,neuropil,info,lights_on,0,0);
temp_data = roifile{i}.Data - artifact(iplane==i,:);
temp_neuropil = roifile{i}.Neuropil - artifact(iplane==i,:);

tex = 4100;
figure
hold on
plot(artifact(iplane==i,tex))
% plot(10*roiline(iplane==i))
plot(diff(temp_data(:,tex-1:tex),[],2))
hold off

%%
figure
iroi = 1;
hold on
% plot(data(iroi,:)-5e3*affected(iroi,:))2
plot(neuropil(iroi,:)-5e3*affected(iroi,:))
% plot(7e3*affected(iroi,:))
hold off

%%
iroi = 1;
figure
hold on
plot(mskcdf(iroi,:))
plot([512*(iplane(iroi)-1)+roiline(iroi) 512*(iplane(iroi)-1)+roiline(iroi)],[0 1])
hold off

%%
figure
hold on
plot(affected(iroi,:))
plot(0.1+affected_r(iroi,:))
hold off

%%
dm = find(diff(affected(iroi,:)));
dr = find(diff(affected_r(iroi,:)));
%%
figure
n = 100;
x = affected(1:n:end);
y = affected_r(1:n:end);
noise = 0.01;
scatter(x+noise*randn(size(y)),y+noise*randn(size(y)))

%%
figure
for iroi=1:30
    subplot(3,10,iroi)
    noise = 0.01;
    x = affected(iroi,1:100:end);
    y = neuropil(iroi,1:100:end);
    scatter(x+noise*randn(size(x)),y,'.')
end

%%
tex = 54;
figure
hold on
plot(neuropil(:,tex))
plot(5e3*affected(:,tex))
hold off

% %
% stimbase = '/home/mossing/modulation/visual_stim/190712/M0148/';
% foldbase = '/home/mossing/modulation/matfiles/190712/M0148/ot/';
% filebase = 'M0148_110_003';
%
% stimbase = '/home/mossing/modulation/visual_stim/190715/M0208/';
% foldbase = '/home/mossing/modulation/matfiles/190715/M0208/ot/';
% filebase = 'M0208_130_004';

% stimbase = '/home/mossing/modulation/visual_stim/190807/M0153/';
% foldbase = '/home/mossing/modulation/matfiles/190807/M0153/ot/';
% sbxbase = '/home/mossing/modulation/2P/190807/M0153/';
% filebase = 'M0153_090_004';

%%

stimbase = '/home/mossing/modulation/visual_stim/191017/M0276/';
foldbase = '/home/mossing/modulation/matfiles/191017/M0276/ot/';
sbxbase = '/home/mossing/modulation/2P/191017/M0276/';
filebase = 'M0276_120_003';

load(sprintf('%s%s.mat',stimbase,filebase),'result');
load(sprintf('%s%s.mat',foldbase,filebase),'info');
for i=1:4
    filenames{i} = sprintf('%s%s_ot_%03d.rois',foldbase,filebase,i-1);
end
%%
while min(diff(info.frame))<0
    ind = find(diff(info.frame)<0,1);
    info.frame(ind+1:end) = info.frame(ind+1:end) + 65536;
end
%%
lights_on = result.gratingInfo.lightsOn;
frm = floor(info.frame(2:4:end-1)/4);
onwhich = rem(info.frame(2:4:end-1)-1,4)+1;
nbefore = 30;
nafter = 50;
avg_trace = zeros(2,nbefore+nafter+1);
for i=1:4
    roifile = load(filenames{i},'-mat');
    for t=1:numel(frm)
        for ll=1:2
            if lights_on(t)==ll-1
                if onwhich(t)~=i
                    avg_trace(ll,:) = avg_trace(ll,:) + nanmean(roifile.corrected(:,frm(t)-nbefore:frm(t)+nafter))/numel(frm);
                end
            end
        end
    end
end
%%
figure
plot(avg_trace')

%%
lights_on = result.gratingInfo.lightsOn>0;
frm = floor(info.frame(2:4:end-1)/4);
%%
onwhich = rem(info.frame(2:4:end-1)-1,4)+1;
avg_trace = cell(4,1); 
for i=1:4
    avg_trace{i} = zeros(2,numel(frm)/2,81);
    roifile = load(filenames{i},'-mat');
    ctr = [0 0];
    for t=1:numel(frm)
        for ll=1:2
            if lights_on(t)==ll-1
                ctr(ll) = ctr(ll)+1;
                avg_trace{i}(ll,ctr(ll),:) = nanmean(roifile.corrected(:,frm(t)-30:frm(t)+50))/numel(frm);
            end
        end
    end
end
%%
figure
hold on
for i=1:4
    subplot(1,4,i)
    plot(squeeze(avg_trace{1}(2,onwhich(lights_on)==i,32)))
end
plot(onwhich(lights_on)==1)
subplot(1,2,1)
for i=1:4
    subplot(1,4,i)
    hist(avg_trace{1}(2,onwhich(lights_on)==i,32),100)
end
subplot(1,2,2)
hist(avg_trace{1}(2,:,33),100)
hold off
%%
nskip = 1;
nbefore = 100; 
nafter = 200;
nboundary = 100;
winaroundtrig = 20;
naround = 3;
nlines = info.sz(1);
nplanes = info.otparam(end);
noffset = nlines*(naround-1)/2;

     
ts_on = 4*(find(lights_on)-1)+1+nskip;
signal = zeros(numel(ts_on),naround*nlines);
for i=1:numel(ts_on)
    im1 = sbxreadpacked(sprintf('%s%s',sbxbase,filebase),info.frame(ts_on(i))-1,naround);
    im2 = sbxreadpacked(sprintf('%s%s',sbxbase,filebase),info.frame(ts_on(i))-1-nplanes,naround);
    signal(i,:) = reshape(mean(im1(:,nboundary+1:end-nboundary,:)-im2(:,nboundary+1:end-nboundary,:),2),1,[]);
end
%%
for i=1:numel(ts_on)
    dif = diff(signal(i,(noffset+info.line(ts_on(i))-winaroundtrig):(noffset+info.line(ts_on(i))+winaroundtrig)));
    dif = dif(1:end-1) + dif(2:end);
    [~,maxind] = max(abs(dif));
    center_on(i) = noffset+info.line(ts_on(i))+maxind-winaroundtrig+1;
    trigaligned(i,:) = signal(i,center_on(i)-nbefore:center_on(i)+nafter);
end
%%
ts_off = nplanes*(find(lights_on)-1)+1+nskip+3;
signal = zeros(numel(ts_off),naround*nlines);
for i=1:numel(ts_off)
    im1 = sbxreadpacked(sprintf('%s%s',sbxbase,filebase),info.frame(ts_off(i))-1,naround);
    im2 = sbxreadpacked(sprintf('%s%s',sbxbase,filebase),info.frame(ts_off(i))-1-nplanes,naround);
    signal(i,:) = reshape(mean(im1(:,nboundary+1:end-nboundary,:)-im2(:,nboundary+1:end-nboundary,:),2),1,[]);
end
%%
for i=1:numel(ts_off)
    dif = diff(signal(i,(noffset+info.line(ts_off(i))-winaroundtrig):(noffset+info.line(ts_off(i))+winaroundtrig)));
    dif = dif(1:end-1) + dif(2:end);
    [~,maxind] = max(abs(dif));
    center_off(i) = noffset+info.line(ts_off(i))+maxind-winaroundtrig+1;
end
%%
t0 = nbefore;
expfun = @(t,t0,a,b,c,d)(t>=t0).*(a*exp(-(t-t0)/b)+d)+c;
expfun_ = @(t,t0,x) expfun(t,t0,x(1),x(2),x(3),x(4));
t = [1:size(trigaligned,2)];
lsqfun = @(model)sum(abs(mean(trigaligned)-model).^1);
costfun = @(x)lsqfun(expfun_(t,t0,x));
c = mean(trigaligned(:,1));
d = mean(trigaligned(:,end))-c;
a = 5e3;
b = 30;
x0 = [a b c d];
xstar = fminunc(costfun,x0);
cstar = costfun(xstar);
% %%
% figure
% hold on
% plot(mean(trigaligned))
% plot(expfun_(t,t0s(imin),xstar(imin,:)))
% % plot(expfun_(t,t0,x0))
% hold off
% %%
% %%
% figure
% hold on
% % plot(mean(trigaligned))
% artifact = repmat(expfun_(t,t0s(imin),xstar(imin,:)),size(trigaligned,1),1);
% plot(trigaligned'-artifact')
% % plot(expfun_(t,t0,x0))
% hold off
% %%

%%

artifact = uint16(zeros(nlines,max(info.frame)));
for i=1:numel(ts_on)
    frame_on = info.frame(ts_on(i));
    frame_off = info.frame(ts_off(i));
    t = 1:nlines*(frame_off-frame_on)+center_off(i)-center_on(i);
    artifact(nlines*(frame_on-1)+center_on(i)-noffset:nlines*(frame_off-1)+center_off(i)-noffset-1) = expfun_(t,0,xstar);
end

artifact = artifact';

%%
%%
figure
hold on
% plot(mean(trigaligned))
itrial = 4;
artifact = repmat(expfun_(t,t0s(imin),xstar(imin,:)),size(trigaligned,1),1);
plot(trigaligned(:,itrial)-artifact(:,itrial))
% plot(expfun_(t,t0,x0))
hold off
%%
i = 400;
npix = size(im1,2);
artifact = zeros(nlines*naround,npix);
t = [1:size(artifact,1)]';
artifact = uint16(repmat(expfun_(t,center_on(i),xstar),1,size(artifact,2)));
artifact = reshape(artifact,[nlines,naround,npix]);
artifact = permute(artifact,[1 3 2]);
iframe = 1;
imshow(im1(:,:,iframe)-uint16(artifact(:,:,iframe)))
