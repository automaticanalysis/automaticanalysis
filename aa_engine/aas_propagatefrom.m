% aa component used by parallel engine
%
%  The directory structure is cached on each machine, so when a file
%  is written it doesn't appear for 8-60 seconds on other machines
%  This made flag passing slow. This uses scp to copy a file to the
%  destination machine so it is available immediately
%
% Rhodri Cusack MRC CBU Aug 2007

function [destfn]=aas_propagatefrom(destmachine, fn, suffix)
if (exist('suffix','var'))
    [pth nme ext]=fileparts(fn);
    destfn=fullfile(pth,[nme suffix ext]);
else
    destfn=fn;
end;

% wait until size stops changing for 1s (indicative of end of file write)
previoussize=-1;
while(1)
    cmd=['ssh -x ' strtrim(destmachine) ' "ls -l ' fn '"'];
    [s w]=unix(cmd);
    if (s) 
        break;
    end;
    [tmp rem]=strtok(w);
    for i=1:4
        [tmp rem]=strtok(rem);
    end;
    sizenow=str2num(tmp);
    if (sizenow==previoussize)
        break;
    end;
    previoussize=sizenow;
    pause(1);
end;

cmd=['scp ' destmachine ':' fn ' '  destfn];
[s w]=unix(cmd);

