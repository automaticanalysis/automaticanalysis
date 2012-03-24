% AA module - copy image orientation
% Copies image orientation from the EPI stream to another stream
% Rhodri Cusack MRC CBU Cambridge Dec 2010

function [aap,resp]=aamod_mapstreams(aap,task,i,j)

resp='';

switch task
    
    case 'report'
    case 'doit'
        
        % Reorient each stream in turn
        inputstreams=aap.tasklist.currenttask.inputstreams.stream;
        outputstreams=aap.tasklist.currenttask.outputstreams.stream;
        if (length(inputstreams)~=length(outputstreams))
            aas_log(aap,true,'Must have same number of input and output streams - check aamod_mapstreams.xml');
        end;
        for streamind=1:length(inputstreams)
            fns= aas_getfiles_bystream(aap,i,j,inputstreams{streamind});
            aap=aas_desc_outputs(aap,i,j,outputstreams{streamind},fns);
        end;
        
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;

