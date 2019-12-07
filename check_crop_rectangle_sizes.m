%%
[foldname,filenames1ch,filenames2ch] = read_exptlist('exptlist.txt');
%%
suite2p_fold = '/home/mossing/data/suite2P/results/';
sbx_fold = '/home/mossing/data/matfiles/'; %2P/';
%%
szs = {};
sz2p = {};
D = {};
for iexpt=1:numel(foldname)
    disp(foldname{iexpt})
    s = strsplit(foldname{iexpt},'/');
    animalid = s{2};
    date = s{1};
    fold_of_interest = [suite2p_fold animalid '/' date '/'];
    d = dirnames(fold_of_interest,fold_of_interest);
    D{iexpt} = d;
    szs{iexpt} = zeros(numel(d),2);
    for ianalysis=1:numel(d)
        disp(['suite2p: ' d{ianalysis}])
        load([d{ianalysis} '/suite2p/plane0/Fall.mat'],'ops')
        szs{iexpt}(ianalysis,:) = size(ops.meanImg);
    end
    disp(['2p: '])
    fold_of_interest = [sbx_fold date '/' animalid '/ot/'];% '/ot/'];
    d = dirnames([fold_of_interest '/M*.mat'],fold_of_interest);
    sz2p{iexpt} = zeros(numel(d),2);
    for imat=1:numel(d)
        load(d{imat},'info')
        try
            sz2p{iexpt}(imat,:) = [info.rect(2)-info.rect(1)+1 info.rect(4)-info.rect(3)+1];
        catch
            disp('')
        end
    end
end