% Automatic analysis - check version
% specify min and max version numbers
%
% Modified to handle _devel on end of version 20/2/2006

function aas_requiresversion(aap)

minver=aap.options.aa_minver;
maxver=aap.options.aa_maxver;

currver=aas_getversion();

if (currver>=minver & currver<=maxver) 
    aas_log(aap,0,sprintf('Current aa version %3.2f suitable for this user script\n',currver));
else
    aas_log(aap,1,sprintf('aa version %f outside of range %f to %f suitable for user script\n',currver,minver,maxver));
end;