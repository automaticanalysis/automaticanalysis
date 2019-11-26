function M = meminfo
%
% return real and virtual memory parameters
% all values returned are in KB
%
% CHANGE HISTORY
%
% [TA] -- added system, limits, and process sections
% [MSJ] -- added OS X version
%

if ismac
	
	% Mac unix (Darwin) and Linux do memory differently -- Linux has
	% a /proc directory and the info can be extracted from there
	% whereas in OS X you have to parse it out of sysctl and ps.
	
	% system

	[ status, result ] = system('sysctl -n hw.memsize');
	M.ResTotal = str2double(result) / 1024; % convert to KB
	
	% this assumes a page size of 4096 bytes:
	
	[ status, result ]= system('vm_stat | grep free');
	spaces = strfind(result,' ');
	M.ResFree = str2num(result(spaces(end):end))*4096 / 1024;
	
	% need to parse SwapTotal and SwapFree from sysctl result
	% -e should force equal sign delimited entries
	
	[ status, result ] = system('sysctl -e vm.swapusage');
		
	result = erase(result,' ');
	
	iStart = strfind(result,'total=');
	temp = result(iStart+length('total='):end);
	M.SwapTotal = sscanf(temp,'%f');
	
	iStart = strfind(result,'free=');
	temp = result(iStart+length('free='):end);
	M.SwapFree = sscanf(temp,'%f');
	
	% vm.swapusage returns results in MB -- convert to KB
	
	M.SwapTotal = M.SwapTotal * 1024;
	M.SwapFree = M.SwapFree * 1024;
	
	% limits

    try
        [ status, result ] = system('ulimit -m');
        result = deblank(result);
        if strcmp(result,'unlimited')
            M.ResLimit = Inf;
        else
            M.ResLimit = str2double(result);
        end
    catch
        fprintf('ulimit not available, defaulting to inf\n');
        M.ResLimit = Inf;
    end

    try
        [ status, result ] = system('ulimit -v');
        result = deblank(result);
        if strcmp(result,'unlimited')
            M.VirtLimit = Inf;
        else
            M.VirtLimit = str2double(result);
        end
    catch
        M.VirtLimit = Inf;
    end
    
	% process
	
	% assuming ps returns results in KB

	[ status, result ] = system(char(sprintf('ps -o rss= %d',feature('getpid'))));
	M.ResUsedBy = str2double(result);
	
	% there's no easy way to get VmSwap for OS X so this currently
	% sets it to zero. The variable is currently unused by any
	% caller and might never be, at least when runnning OS X.

	M.SwapUsedBy = 0;
	
	[ status, result ] = system(char(sprintf('ps -o vsz= %d',feature('getpid'))));
	M.VirtUsedBy = str2double(result);
		
else

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

    % ulimit is bash specific so this will fail on csh / tcsh
    [exception, out] = aas_shell('ulimit -m', true, false);
    if exception
        fprintf('ulimit not available, defaulting to inf\n');
        dat = Inf;
    else
        out = deblank(out);
        if strcmp(out,'unlimited'), dat = Inf;
        else dat = str2double(out); end
    end
	M.ResLimit = dat;

    [exception, out] = aas_shell('ulimit -v',true, false);
    if exception
        dat = Inf;
    else
        out = deblank(out);
        if strcmp(out,'unlimited'), dat = Inf;
        else dat = str2double(out); end
    end
	M.VirtLimit = dat;

	% process

	[exception, out] = aas_shell(sprintf('cat /proc/%d/status',feature('getpid')),true);
	out = textscan(out,'%s','delimiter',':'); out = out{1};
	dat = str2double(strtok(out{find(strcmp(out,'VmRSS'))+1}));
	M.ResUsedBy = dat;
	dat = str2double(strtok(out{find(strcmp(out,'VmSwap'))+1}));
	M.SwapUsedBy = dat;
	dat = str2double(strtok(out{find(strcmp(out,'VmSize'))+1}));
	M.VirtUsedBy = dat;


end








end

