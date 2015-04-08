function dep = dep_read(fname)
if ~nargin
    fname = 'aap_prov.trp';
end
if ~exist(fname,'file'), error('Provanence %s not found!',fname); end
fid = fopen(fname,'r');
% dep = {};

while ~feof(fid)
    line = fgetl(fid);
    [src, res] = strtok_ptrn(line(2:end),'> <');
    [str, res] = strtok_ptrn(res(4:end),'> <');
    if any(str=='.')
        [junk, str] = strtok(str,'.');
        str = str(2:end);
    end
    trg = res(4:end-3);
    dep.(trg).(str) = src;
    %     dep(end+1,1:3) = {trg str src};
end

fclose(fid);
end