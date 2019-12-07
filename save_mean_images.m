function save_mean_images(fname)
%%
% figure out which experiments to look at
cd /home/mossing/notebooks/opto1p
[foldname,filenames1ch,filenames2ch] = read_exptlist(fname);

sbx_base = '/media/data/dan/2P/' % this is the location which the mean_images 
% files will be saved to

%%
% lkat = [0 0 1 1 1 1 0]<1;
% foldname = foldname(lkat);
% filenames1ch = filenames1ch(lkat);
% filenames2ch = filenames2ch(lkat);

%%
% figure out suite2p filenames
result_base1ch = '/home/mossing/data1/suite2P/results/';
basefolds1ch = gen_suite2p_filenames(foldname,filenames1ch,result_base1ch);
result_base2ch = '/home/mossing/data_ssd/suite2P/results/';
basefolds2ch = gen_suite2p_filenames(foldname,filenames2ch,result_base2ch);
% result_base = '/home/mossing/data/suite2P/results/';
% [basefolds,animalid,date] = deal(cell(size(foldname)));
for isesh=1:numel(foldname)
    fileparts = strsplit(foldname{isesh},'/');
    animalid{isesh} = fileparts{2};
    date{isesh} = fileparts{1};
end
%     s = cell(size(filenames1ch{isesh}));
%     for isubfold = 1:numel(filenames1ch{isesh})
%         s{isubfold} = num2str(filenames1ch{isesh}(isubfold));
%     end
%     subfold = strjoin(s,'_');
%     basefolds{isesh} = strjoin({result_base, animalid{isesh}, date{isesh}, subfold, 'suite2p'},'/');
% end

%%
% % load mean image from 1-channel analysis
% for i=1:numel(basefolds1ch)
%     
% end

%%
% save mean images into suite2p file structure
% NEED TO ADD DEFORMATION
for isesh=1:numel(basefolds2ch)
    basefold = basefolds2ch{isesh};
    old_fold = pwd;
    cd(basefold)
    disp(basefold)
    for iplane=0:3
        planefoldname = ['plane' num2str(iplane)];
        fold1ch = [basefolds1ch{isesh} '/' planefoldname];
        im = load_suite2p_mean_image(fold1ch);
        %im_compare: 1-channel data
        save_suite2p_mean_images(planefoldname,im)
        close('all')
%         cd(basefold)
    end
    cd(old_fold)
end

%%
% copy mean images into sbx file structure
for isesh=1:numel(basefolds2ch)
    targetfold = strjoin({sbx_base, date{isesh}, animalid{isesh}},'/');
    if ~exist(targetfold)
        mkdir(targetfold)
    end
    for iplane=0:3
        planefoldname = ['plane' num2str(iplane)];
        sourcefile = [basefolds2ch{isesh} '/' planefoldname '/mean_images.mat'];
        targetfile = sprintf([targetfold '/mean_images_ot_%03d.mat'],iplane);
        copyfile(sourcefile,targetfile);
    end
end

function basefolds = gen_suite2p_filenames(foldname,filenames,result_base)
% figure out suite2p filenames of 1ch data

[basefolds,animalid,date] = deal(cell(size(foldname)));
for isesh=1:numel(foldname)
    fileparts = strsplit(foldname{isesh},'/');
    animalid{isesh} = fileparts{2};
    date{isesh} = fileparts{1};
    s = cell(size(filenames{isesh}));
    for isubfold = 1:numel(filenames{isesh})
        s{isubfold} = num2str(filenames{isesh}(isubfold));
    end
    subfold = strjoin(s,'_');
    basefolds{isesh} = strjoin({result_base, animalid{isesh}, date{isesh}, subfold, 'suite2p'},'/');
end
