% aa module tempate for sessions 

function [aap,resp]=aamod_session(aap,task,subj,sess)
resp='';

switch task
    case 'report'
        localpath = aas_getpath_bydomain(aap,aap.tasklist.currenttask.domain,[subj,sess]);
        
        fdiag = dir(fullfile(localpath,'diagnostic_*.jpg'));
        if isempty(fdiag)
            streams=aas_getstreams(aap,'output');
            for streamind=1:length(streams)
                % obtain output
                outputfnames = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,[subj sess],streams{streamind},'output');
                
                % perform diagnostics
                do_diag(outputfnames);
            end
            fdiag = dir(fullfile(localpath,'diagnostic_*.jpg'));
        end
        
        for d = 1:numel(fdiag)
            aap = aas_report_add(aap,subj,'<table><tr><td>');
            imgpath = fullfile(localpath,fdiag(d).name);
            aap=aas_report_addimage(aap,subj,imgpath);
            aap = aas_report_add(aap,subj,'</td></tr></table>');
        end
    case 'doit'
        
        streams=aas_getstreams(aap,'input');
        
        for streamind=1:length(streams)
            % obtain input filename(s)
            inputfnames = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,[subj sess],streams{streamind});
            pfx = 'o'; % prefix for outputs
            outputfnames = spm_file(inputfnames,'prefix',pfx);
            
            % obtain parameter
            par   = aap.tasklist.currenttask.settings.parameter;
            
            % do stuff
            do_stuff(inputfnames, outputfnames, par);
            
            % Describe outputs
            aap=aas_desc_outputs(aap,aap.tasklist.currenttask.domain,[subj sess],streams{streamind},outputfnames);
            
        end;
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end