%%
cd /home/mossing/notebooks/opto1p
[foldname,filenames] = read_exptlist(fname);

%%
result_base = '/home/mossing/data/suite2P/results/';
sbx_base = '/home/mossing/data/2P/';
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

%%
for isesh=1:numel(basefolds)
    basefold = basefolds{isesh};
    old_fold = pwd;
    cd(basefold)
    for iplane=0:3
        planefoldname = ['plane' num2str(iplane)];
        save_suite2p_mean_images(planefoldname)
    end
    cd(old_fold)
end

%%
for isesh=1:numel(basefolds)
    targetfold = strjoin({sbx_base, date{isesh}, animalid{isesh}},'/');
    if ~exist(targetfold)
        mkdir(targetfold)
    end
    for iplane=0:3
        planefoldname = ['plane' num2str(iplane)];
        sourcefile = [basefolds{isesh} '/' planefoldname '/mean_images.mat'];
        targetfile = sprintf([targetfold '/mean_images_ot_%03d.mat'],iplane);
        copyfile(sourcefile,targetfile);
    end
end