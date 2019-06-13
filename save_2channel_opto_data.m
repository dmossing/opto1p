% %%
% addpath(genpath('/home/mossing/Documents/code/s2p_current'));
% addpath(genpath('/home/mossing/Documents/code/adesnal'));
% addpath(genpath('/home/mossing/Documents/code/sbatch_scripts'));
% addpath(genpath('/home/mossing/Documents/code/downloads/saveastiff'));
% addpath(genpath('~/Documents/code/OASIS_matlab'));addpath(genpath('~/downloads/sort_nat'));
% addpath(genpath('/home/mossing/Documents/code/downloads/EvansCode/'))
% addpath(genpath('/home/mossing/Documents/code/downloads/npy-matlab-master'))

%%

% FIRST RUN save_and_transfer_crop(foldname)

% foldname = '/home/mossing/data/2P/181030/M9826/';
% foldname = '/home/mossing/data/2P/190102/M10130/';
% foldname = '/home/mossing/modulation/2P/181213/M10345/';
% foldname = '/home/mossing/modulation/2P/181214/M10130/';
% foldname = '/home/mossing/data/2P/181117/M10039/';
% foldname = '/home/mossing/data/2P/181120/M010039/';
% foldname = '/home/mossing/data/2P/181205/M10130/';

% foldname = '/media/mossing/data_ssd/data/2P/green_only/190221/M9835/';
% foldname = '/media/mossing/data_ssd/data/2P/190225/M10344/';
foldname = {};
filenames = {};
k = 1;

foldname{k} = '190301/M9835/';
filenames{k} = [1,888];
k = k+1;

foldname{k} = '190307/M9835/';
filenames{k} = [1,999];
k = k+1;

foldname{k} = '190318/M10338/';
filenames{k} = [2,999];
k = k+1;

foldname{k} = '190320/M10365/';
filenames{k} = [1,999];
k = k+1;

foldname{k} = '190410/M10368/';
filenames{k} = [2];
k = k+1;

foldname{k} = '190411/M0002/';
filenames{k} = [4];
k = k+1;

foldname{k} = '190501/M0094/';
filenames{k} = [4];
k = k+1;

%%

% data_foldbase = '/media/mossing/backup_0/data/2P/';
% result_foldbase = '/home/mossing/modulation/visual_stim/';
data_foldbase = '/media/greg/modulation/mossing/2P/';
result_foldbase = '/media/greg/modulation/mossing/visual_stim/';

targetfold = '/media/mossing/data/suite2P/raw/';

opts.chunksize = 1000;
opts.green_only = 0;


for k=1:length(foldname)
    thisfoldname = foldname{k};
    d = dir([data_foldbase thisfoldname '/M*.mat']);
    fnames = {d(:).name};
    for i=1:numel(d)
        s = strsplit(fnames{i}(1:end-4),'_');
        exptno = str2num(s{end});
        if ismember(exptno,filenames{k})
            fileparts = strsplit(thisfoldname,'/');
            animalid = fileparts{2};
            dstr = fileparts{1};
            subfold = num2str(exptno);
            opts.targetfold = [targetfold animalid '/' dstr '/' subfold '/'];
            sbx_to_cropped_tiffs([data_foldbase thisfoldname '/' fnames{i}(1:end-4)],opts); % 1 for green only; otherwise 0
%             move_suite2p_tiffs([data_foldbase thisfoldname '/' fnames{i}(1:end-4)],targetfold,'2P');
        end
    end
    
end