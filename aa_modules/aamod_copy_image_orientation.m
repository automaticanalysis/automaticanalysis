% AA module - copy image orientation
% Copies image orientation from the first stream to the other streams
% Rhodri Cusack MRC CBU Cambridge Dec 2010

function [aap,resp]=aamod_copy_image_orientation(aap,task,i,j)

resp='';

switch task
    
    case 'report'
    case 'doit'
        
        
        streams=aap.tasklist.currenttask.inputstreams.stream;

        % get first stream
        fns_epi= aas_getfiles_bystream(aap,i,j,streams{1});
        
        % Reorient each stream in turn
        for streamind=2:length(streams)  
            fns_out=aas_getfiles_bystream(aap,i,j,streams{streamind});
            if (size(fns_epi,1)~=size(fns_out,1))
                aas_log(aap,true,sprintf('Stream epi and stream %s must contain exactly the same number of files, but got %d and %d.',streams{streamind},size(fns_epi,2),size(fns_out,2)));
            end;
            for fileind=1:size(fns_epi,1)
                spm_get_space(fns_out(fileind,:),spm_get_space(fns_epi(fileind,:)));
            end;
            aap=aas_desc_outputs(aap,i,j,streams{streamind},fns_out);
        end;
        

case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;

