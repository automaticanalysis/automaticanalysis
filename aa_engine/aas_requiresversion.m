% Automatic analysis - check version
% specify min and max version numbers
%
% Modified to handle _devel on end of version 20/2/2006

function aas_requiresversion(aap)

global aa
currver = sscanf(aa.Version,'%d.%d.%d')'*[100^2 100 1]';
minver  = sscanf(aap.options.aa_minver,'%d.%d.%d')'*[100^2 100 1]';
maxver  = sscanf(aap.options.aa_maxver,'%d.%d.%d')'*[100^2 100 1]';

if (currver>=minver && currver<=maxver) 
    aas_log(aap,0,sprintf('Current aa version %3.2f suitable for this user script\n',currver));
else
    aas_log(aap,1,sprintf('aa version %f outside of range %f to %f suitable for user script\n',currver,minver,maxver));
end;