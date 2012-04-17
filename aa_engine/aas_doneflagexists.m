%%==============================
% Check for done flag
function [resp]=aas_doneflagexists(aap,doneflag)
global aaworker
switch(aap.directory_conventions.remotefilesystem)
    case 's3'
        [aap attr]=sdb_get_attributes(aap,aaworker.doneflagtablename,doneflag);
        resp=~isempty(attr);
%         if (resp)
%             aas_log(aap,false,sprintf('Found done flag %s in table %s',doneflag, aaworker.doneflagtablename)); 
%         else
%             aas_log(aap,false,sprintf('No done flag %s in table %s',doneflag, aaworker.doneflagtablename)); 
%         end;
    otherwise
        resp=exist(doneflag,'file');
end;
