function M = meminfo
% system
[junk, out] = aas_shell('cat /proc/meminfo',true);
out = textscan(out,'%s','delimiter',':'); out = out{1};
dat = str2double(strtok(out{find(strcmp(out,'MemTotal'))+1}));
M.ResTotal = dat;
dat = str2double(strtok(out{find(strcmp(out,'MemFree'))+1}));
M.ResFree = dat;
dat = str2double(strtok(out{find(strcmp(out,'SwapTotal'))+1}));
M.SwapTotal = dat;
dat = str2double(strtok(out{find(strcmp(out,'SwapFree'))+1}));
M.SwapFree = dat;

% limits
[junk, out] = aas_shell('ulimit -m',true); out = deblank(out);
if strcmp(out,'unlimited'), dat = Inf;
else dat = str2double(out); end
M.ResLimit = dat;
[junk, out] = aas_shell('ulimit -v',true); out = deblank(out);
if strcmp(out,'unlimited'), dat = Inf;
else dat = str2double(out); end
M.VirtLimit = dat;

% process
[junk, out] = aas_shell(sprintf('cat /proc/%d/status',feature('getpid')),true);
out = textscan(out,'%s','delimiter',':'); out = out{1};
dat = str2double(strtok(out{find(strcmp(out,'VmRSS'))+1}));
M.ResUsedBy = dat;
dat = str2double(strtok(out{find(strcmp(out,'VmSwap'))+1}));
M.SwapUsedBy = dat;
dat = str2double(strtok(out{find(strcmp(out,'VmSize'))+1}));
M.VirtUsedBy = dat;
end

