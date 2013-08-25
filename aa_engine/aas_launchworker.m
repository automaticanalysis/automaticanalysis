function [s]=aas_launchworker(aap,workernumber,specialrequirements)
if (~exist('specialrequirements','var'))
    specialrequirements={};
end;

pth=aaworker_getparmpath(aap,workernumber);
[s rem]=unix('hostname');
% Just take the last line in case there are warnings (occasionally get
% Xauth error on our system)
rem=strtrim(rem);
while(1)
    [hostname rem]=strtok(rem,10);
    if (isempty(rem))
        break;
    end;
end;
hostname=strtrim(hostname);
rmdir(pth,'s');
aas_log(aap,0,sprintf('PARALLEL launching worker %d at %s',workernumber,datestr(now,14)));
parmpath=aaworker_getparmpath(aap,workernumber);
savepath(fullfile(parmpath,'linemanagerpaths.m'));
specialreqstr='';
for j=1:length(specialrequirements)
    specialreqstr=[specialreqstr ' ' specialrequirements{j}];
end;
cmd=sprintf('p=linemanagerpaths; path(p); aaworker_prepare(%d,[%s])',workernumber,sprintf('%d ',double(hostname)));
cmdall=['cd  ' parmpath '; spm title:aaworker' num2str(workernumber) specialreqstr ' workerdesktop fmri "' cmd '"'];
[s w]=aas_shell(cmdall);
if (s)
    aas_log(aap,false,sprintf('AAPARALLEL warning: worker %d said:\n%s',workernumber,w));
end;
