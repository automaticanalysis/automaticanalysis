% AA module - first level statistics 2: convolve model (spm_fmri_spm_ui)
%
% Split version of aamod_firstlevel_model - useful for inserting custom model
% options in between stages. Replaces aamod_firstlevel_model as a whole. Stages:
%
% aamod_firstlevel_model_1_config - add events, set up basic SPM struct
% aamod_firstlevel_model_2_convolve - generate convolved design matrix
% (spm_fmri_spm_ui)
% aamod_firstlevel_model_3_estimate - fit the model (spm_spm)

function [aap,resp]=aamod_firstlevel_model_2_convolve(aap,task,subj)

resp='';

switch task
    case 'report' % [TA]
        if ~exist(fullfile(aas_getsubjpath(aap,subj),'diagnostic_aamod_firstlevel_model_design.jpg'),'file')
            load(aas_getfiles_bystream(aap,subj,aap.tasklist.currenttask.outputstreams.stream{1}));
            spm_DesRep('DesOrth',SPM.xX);
            saveas(spm_figure('GetWin','Graphics'),fullfile(aas_getsubjpath(aap,subj),'diagnostic_aamod_firstlevel_model_design.jpg'));
            close all;
        end
        fdiag = dir(fullfile(aas_getsubjpath(aap,subj),'diagnostic_*.jpg'));
        for d = 1:numel(fdiag)
            aap = aas_report_add(aap,subj,'<table><tr><td>');
            aap=aas_report_addimage(aap,subj,fullfile(aas_getsubjpath(aap,subj),fdiag(d).name));
            aap = aas_report_add(aap,subj,'</td></tr></table>');
        end

    case 'doit'

        cwd=pwd;
        spmpath = aas_getfiles_bystream(aap,subj,'firstlevel_spm');
        load(spmpath);
        anadir = fileparts(spmpath);
        cd(anadir)
        SPM.swd = anadir;

        %%%%%%%%%%%%%%%%%%%
        %% DESIGN MATRIX %%
        %%%%%%%%%%%%%%%%%%%
        if aap.tasklist.currenttask.settings.firstlevelmasking && aap.tasklist.currenttask.settings.firstlevelmasking < 1
            spm_get_defaults('mask.thresh',aap.tasklist.currenttask.settings.firstlevelmasking);
        end
        SPMdes = spm_fmri_spm_ui(SPM);

        SPMdes.xX.X = double(SPMdes.xX.X);
        
        % DIAGNOSTIC
        subjname = aas_prepare_diagnostic(aap, subj);
        try
            saveas(1, fullfile(aap.acq_details.root, 'diagnostics', ...
                                                [mfilename '__' subjname '.fig']));
        catch
        end

        % now check real covariates and nuisance variables are
        % specified correctly (by using hand-set variables from
        % aamod_firstlevel_model_1_config).
        SPMdes.xX.iG = SPM.xX.iG;
        SPMdes.xX.iC = SPM.xX.iC;

        % Turn off masking if requested
        if ~aap.tasklist.currenttask.settings.firstlevelmasking
            SPMdes.xM=-inf(size(SPMdes.xX.X,1),1);
        end
        
        % Correct epmty model for sphericity check
        if isempty([SPMdes.Sess.U])
            SPMdes.xX.W  = sparse(eye(size(SPMdes.xY.P,1)));
            SPMdes.xVi.V = sparse(eye(size(SPMdes.xY.P,1)));            
        end
        
        %% Describe outputs
        cd (cwd);

        % Describe outputs
        %  firstlevel_spm
        aap=aas_desc_outputs(aap,subj,'firstlevel_spm',fullfile(anadir,'SPM.mat'));

        if aap.tasklist.currenttask.settings.firstlevelmasking
            mask = dir(fullfile(anadir,'mask.*'));
            mask=strcat(repmat([anadir filesep],[numel(mask) 1]),char({mask.name}));
            aap=aas_desc_outputs(aap,subj,'firstlevel_brainmask',mask);            
        end
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
