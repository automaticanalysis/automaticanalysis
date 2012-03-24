% Automatic analysis - check version
% Specify min and max version numbers in 
% aap.options.aa_minver and aap.options.aa_maxver 
% and then call this function in your user script.

function aa_requiresversion(aap)

minver=aap.options.aa_minver;
maxver=aap.options.aa_maxver;

currpth=path;
startverpos=findstr('aa_ver',currpth);
endverpos=findstr(filesep,currpth(startverpos(1):length(currpth)));
currver=str2num(currpth(startverpos(1):endverpos(1)));

if (currver>=minver & currver<=maxver) 
    aas_log(aap,0,sprintf('Current aa version %f suitable for user script\n',currver));
else
    aas_log(aap,1,sprintf('aa version %f outside of range %f to %f suitable for user script\n',currver,minver,maxver));
end;