function [foldname,filenames] = read_exptlist(fname)

%%
fid = fopen(fname,'r');
content = textscan(fid,'%s','delimiter','\n');
content = content{1};
fclose(fid);

%%
lines_per_expt = 3;
foldname = {};
filenames = {};
for i=1:numel(content)
    if mod(i,lines_per_expt)==1
        foldname = [foldname(:); content(i)];
    elseif mod(i,lines_per_expt)==2
        filenames = [filenames(:); {str2num(content{i})}];
    end
end