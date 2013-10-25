% AA module - convert 3d analyze to 4d
% Rhodri Cusack UWO Jun 2011
% rhodri@cusacklab.org

function [aap,resp]=aamod_3dto4d(aap,task,varargin)

resp='';

switch task
    case 'report'
        
    case 'doit'
        if length(aap.tasklist.currenttask.inputstreams.stream)~=length(aap.tasklist.currenttask.outputstreams.stream)
            aas_log(aap,true,'Number of input streams must be the same as number of output streams - check .xml accompanying this module');
        end;
        
        % Get inputs
        for sind=1:length(aap.tasklist.currenttask.inputstreams.stream)
            
            subj_imgs=aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,[varargin{:}],aap.tasklist.currenttask.inputstreams.stream{sind});
            
            % Find common part to filenames and use this as the base of the output
            incommon=subj_imgs==repmat(subj_imgs(1,:),[size(subj_imgs,1) 1]);
            commonpart=find(~all(incommon),1,'first')-1;
            if (commonpart<1)
                outputname='xxx_4d';
            else
                outputname=[subj_imgs(1,1:commonpart) '_4d'];
            end;
            subj_imgs(:,end+1)=' ';
            subj_imgs=subj_imgs';
            subj_imgs=subj_imgs(:)';
            cmd=['fslmerge -t ' outputname ' ' subj_imgs];
            [s w]=aas_runfslcommand(aap,cmd);
            
            if (s)
                aas_log(aap,true,sprintf('Error \n  %s\nwhen executing fsl command\n  %s',w,cmd));
            end;
            outfns=dir([outputname,'*']);
            if (length(outfns)~=1)
                aas_log(aap,true,'Problem identifying output of fsl 3d to 4d');
            end;
            
            [pth nme ext]=fileparts(outputname);
            
            aap=aas_desc_outputs(aap,aap.tasklist.currenttask.domain,[varargin{:}],aap.tasklist.currenttask.outputstreams.stream{sind},fullfile(pth,outfns(1).name));
            
        end;
end;