%%==============================
% Delete a done flag by its path

function aas_delete_doneflag_bypath(aap,doneflag)
global aaworker
switch(aap.directory_conventions.remotefilesystem)
    case 's3'
        fprintf('Deleting done flag %s\n',doneflag);
        attr=sdb_delete_attributes(aap,aaworker.doneflagtablename,doneflag);
        
    otherwise
        if (exist(doneflag,'file'))
            cmd=['rm ' doneflag];
            [s w]=aas_shell(cmd);
            if (s~=0)
                aas_log(aap,1,w);
            end;
        end;
end;
