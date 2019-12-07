function [foldname,filenames1ch,filenames2ch] = read_exptlist(fname)

%%
fid = fopen(fname,'r');
content = textscan(fid,'%s','delimiter','\n');
content = content{1};
fclose(fid);

%%
lines_per_expt = 4;
foldname = {};
filenames1ch = {};
filenames2ch = {};
for i=1:numel(content)
    % check if commented
    if ~isempty(content{i}) && ~strcmp(content{i}(1),'#')
        if mod(i,lines_per_expt)==1
            foldname = [foldname(:); content(i)];
        elseif mod(i,lines_per_expt)==2
            filenames1ch = [filenames1ch(:); {str2num(content{i})}];
        elseif mod(i,lines_per_expt)==3
            filenames2ch = [filenames2ch(:); {str2num(content{i})}];
        end
    end
end
assert(numel(filenames1ch)==numel(filenames2ch));
for i=1:numel(filenames1ch)
    filenames1ch{i} = sort([filenames1ch{i} filenames2ch{i}]);
end