% AA module - reslicing
% [aap,resp]=aamod_reslice(aap,task,i)
% Reslices images to a common space defined by the first image
% i=subject num
% Rhodri Cusack MRC CBU 2004-6 based on original by Matthew Brett
% Modified by Rik Henson 2006-8 to accept reslice "which" option
% 	(plus more defaults can be passed)

function [aap,resp]=aamod_reslice(aap,task,i)

resp='';

switch task
    case 'report'
    case 'doit'
        
        streams=aap.tasklist.currenttask.inputstreams.stream;
        for streamind=1:length(streams)
            for j=aap.acq_details.selected_sessions
                imgs{j} = aas_getimages_bystream(aap,i,j,streams{streamind});
            end;
            % Run the reslicing
            spm_reslice(imgs);
            % Describe outputs
            for j = aap.acq_details.selected_sessions
                rimgs=[];
                for k=1:size(imgs{j},1);
                    [pth nme ext]=fileparts(imgs{j}(k,:));
                    rimgs=strvcat(rimgs,fullfile(pth,['r' nme ext]));
                end;
                sessdir=aas_getsesspath(aap,i,j);
                aas_desc_outputs(aap,i,j,streams{streamind},rimgs);
            end;
        end;
        
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;














