% This is a tempate for a module code processing an MRI session

function [aap,resp]=aamod_reorienttomiddle(aap,task,varargin)
resp='';

switch task
    case 'report'
%         localpath = aas_getpath_bydomain(aap,aap.tasklist.currenttask.domain,[subj,sess]);
%         
%         fdiag = dir(fullfile(localpath,'diagnostic_*.jpg'));
%         if isempty(fdiag)
%             streams=aas_getstreams(aap,'output');
%             for streamind=1:length(streams)
%                 % obtain output
%                 outputfnames = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,[subj sess],streams{streamind},'output');
%                 
%                 % perform diagnostics
%                 do_diag(outputfnames);
%             end
%             fdiag = dir(fullfile(localpath,'diagnostic_*.jpg'));
%         end
%         
%         for d = 1:numel(fdiag)
%             aap = aas_report_add(aap,subj,'<table><tr><td>');
%             imgpath = fullfile(localpath,fdiag(d).name);
%             aap=aas_report_addimage(aap,subj,imgpath);
%             aap = aas_report_add(aap,subj,'</td></tr></table>');
%         end
    case 'doit'
        % init
        indices = cell2mat(varargin);        
        streams=aas_getstreams(aap,'input');
        
        % obtain input filename(s)
        fname = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,indices,streams{1});
        
        % do stuff
        M = eye(4);
        bbcentre = mean(eye(4)*[spm_get_bbox(fname)';1 1],2);
        M(1:3,4) = -bbcentre(1:3);
        
        job.srcfiles = cellstr(spm_select('ExtFPListRec', aas_getpath_bydomain(aap,aap.tasklist.currenttask.domain,indices), basename(fname)));
        job.transform.transM = M;
        job.prefix = '';
        
        spm_run_reorient(job)
        
        % Describe outputs
        aap=aas_desc_outputs(aap,aap.tasklist.currenttask.domain,indices,streams{1},fname);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end