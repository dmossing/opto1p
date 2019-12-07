function trace = compute_rowwise_average_from_raw(filebase,sbxbase)
load(sprintf('%s%s',sbxbase,filebase),'info')

% deal with overflow in trigger frame variable
while find(diff(info.frame)<0,1)
    seam = find(diff(info.frame)<0,1);
    info.frame(seam+1:end) = info.frame(seam+1:end)+65536;
end

fr_no = max(info.frame)-1;
row_no = info.sz(1);
nboundary = 100;
chunksize = 1000;
trace = zeros(row_no,fr_no);

for i=1:chunksize:fr_no
    disp(sprintf('%d/%d',i,fr_no))
    thischunksize = min(fr_no-i,chunksize);
    im = sbxreadpacked(sprintf('%s%s',sbxbase,filebase),i-1,thischunksize);
    trace(:,i:i+thischunksize-1) = mean(im(:,nboundary+1:end-nboundary,:),2);
end