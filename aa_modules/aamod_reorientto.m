function [aap,resp]=aamod_reorientto(aap,task,varargin)
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
        
        % target
        M = spm_get_space(aas_getfiles_bystream_multilevel(aap,aap.tasklist.currenttask.domain,indices,streams{1}));
        
        % do stuff
        sourceimfn = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,indices,streams{2});
        sV = spm_vol(sourceimfn);
        for i = 1:numel(sV)
            spm_get_space(spm_file(sV(i).fname,'number',sV(i).n(1)), M);
        end
        
        % Describe outputs
        aap=aas_desc_outputs(aap,aap.tasklist.currenttask.domain,indices,streams{2},sourceimfn);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end