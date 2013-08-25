% Automatic analysis
% Gets process ID of matlab for a particular worker
% There must be a more elegant way to find the process ID within the
% worker, but I cannot find it
% Rhodri Cusack CBU Cambridge July 2008
% function [pid]=aas_getworkerpid(hostname,workerprocesskey)

function [pid]=aas_getworkerpid(aap,hostname,workerprocesskey)

cmd=['ssh -x ' strtrim(hostname) ' "ps -elf | grep aaworker' num2str(workerprocesskey) ' | grep -v grep"'];
[s rem]=unix(cmd);
if (~s)

    % Make a list of the pids if there is more than one
    pid=[];
    while(length(rem)>0)
        [remline rem]=strtok(rem,10);
        for i=1:4
            [part remline]=strtok(remline);
        end;
        if (length(part)>0)
            pid=[pid str2num(part)];
        end;
    end;
end;