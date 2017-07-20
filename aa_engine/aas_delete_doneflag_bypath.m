%%==============================
% Delete a done flag by its path

function aas_delete_doneflag_bypath(aap,doneflag)
global aaworker
switch(aap.directory_conventions.remotefilesystem)
    case 's3'
        aas_log(aap,false,sprintf('Deleting done flag %s',doneflag));
        attr=sdb_delete_attributes(aap,aaworker.doneflagtablename,doneflag);
        
    otherwise
        if exist(doneflag,'file')
            delete(doneflag)
        end;
end;
