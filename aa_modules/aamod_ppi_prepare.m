% AA module - first level statistics for PPI
% [aap,resp]=aamod_ppi_model(aap,task,subj,sess)
% Limitations:
%   - single VOI
%   - single session
%   - no parametric modulation
% First-level model based on ppi_spm_batch.m 17 by Guillaume Flandin & Darren Gitelman
% Tibor Auer MRC CBSU Jul 2014

function [aap,resp]=aamod_ppi_prepare(aap,task,subj,sess)

resp='';

switch task
    case 'domain'
        resp='session';   % this module needs to be run once per session case 'report'
        
    case 'doit'
        
        %% Init
        % Directory
        sess_dir = aas_getsesspath(aap,subj,sess);
        
        % Settings
        VOIs = aas_getfiles_bystream(aap,subj,sess,'vois');
        load(deblank(VOIs(aap.tasklist.currenttask.settings.voi,:)));
        
        %% Prepare
        cd(sess_dir);
        spm_jobman('initcfg');
        fSPM = aas_getfiles_bystream(aap, subj,'firstlevel_spm');
        load(fSPM);
        
        % PPI
        conname = '';
        ppicon = zeros(0,3);
        Ic = aap.tasklist.currenttask.settings.contrasts;
        iAct = find(Ic>0); iDeact = find(Ic<0);
        for i = 1:numel(iAct)
            conname = [conname '+' SPM.Sess(sess).U(iAct(i)).name{1}];
            ppicon = vertcat(ppicon,[iAct(i) 1 1]);
        end
        for i = 1:numel(iDeact)
            conname = [conname '-' SPM.Sess(sess).U(iDeact(i)).name{1}];
            ppicon = vertcat(ppicon,[iDeact(i) 1 -1]);
        end
        ppiname = [strtok(xY.name,'-0') 'x' conname ];
        
        spm_peb_ppi(SPM,'ppi',xY,ppicon,ppiname,0);
        fPPI = spm_select('FPList',sess_dir,'^PPI_.*\.mat');
        
        %% Describe outputs
        %  firstlevel_spm
        aap=aas_desc_outputs(aap,subj,sess,'ppi',fPPI);

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;