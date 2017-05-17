function [aap,resp]=aamod_split(aap,task,subj,sess)
resp='';

switch task
    case 'report'
        
    case 'doit'
        pfx = 'c'; % prefix for outputs
        
        stream=aas_getstreams(aap,'input'); stream = stream{1};
        inputfnames = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,[subj sess],stream);
        
        % split
        Vo = spm_file_split(spm_vol(inputfnames), aas_getsesspath(aap,subj,sess));
        
        % remove volume(s) of no interest
        toRemove = Vo;
        toRemove(aap.tasklist.currenttask.settings.start:min(aap.tasklist.currenttask.settings.stop,numel(Vo))) = [];
        delete(toRemove.fname);
        Vo = Vo(aap.tasklist.currenttask.settings.start:min(aap.tasklist.currenttask.settings.stop,numel(Vo)));
        
        % obtain parameter
        if aap.tasklist.currenttask.settings.NIFTI4D
            finalepis = Vo(1).fname;
            ind = find(finalepis=='_');
            finalepis = spm_file(finalepis(1:ind(end)-1),'prefix',pfx,'ext','nii');
			spm_file_merge(Vo,finalepis,0);
            delete(Vo.fname);
        else
            finalepis = {};
            for v = Vo'
                finalepis{end+1} = spm_file(v.fname,'prefix',pfx);
                movefile(v.fname,finalepis{end});
            end
            finalepis = char(finalepis);
        end        
        
        % Describe outputs
        aap=aas_desc_outputs(aap,aap.tasklist.currenttask.domain,[subj sess],stream,finalepis);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,true,sprintf('Unknown task %s',task));
end