function mem = meminfo
fid = fopen('/proc/meminfo');
dat = textscan(fid,'%s %d kB');
fclose(fid);
for f = 1:numel(dat{1})
    field = dat{1}{f}(1:end-1);
    field = strrep(field,'(','_');
    field = strrep(field,')','');
    mem.(field) = dat{2}(f);
end