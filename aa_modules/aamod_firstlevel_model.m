% AA module - first level statistics
% **********************************************************************
% You should no longer need to change this module - you may just
% modify the .xml or model in your user script
% **********************************************************************
% Based on original by FIL London and Adam Hampshire MRC CBU Cambridge Feb 2006
% Modified for aa by Rhodri Cusack MRC CBU Mar 2006-Aug 2007
% Thanks to Rik Henson for various suggestions (modified) [AVG & TA]

function [aap,resp]=aamod_firstlevel_model(aap,task,subj)

resp='';

switch task
    case 'report' % [TA]
        if ~exist(fullfile(aas_getsubjpath(aap,subj),['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_design.jpg']),'file')
            load(aas_getfiles_bystream(aap,subj,aap.tasklist.currenttask.outputstreams.stream{1}));
            spm_DesRep('DesOrth',SPM.xX);
            saveas(spm_figure('GetWin','Graphics'),fullfile(aas_getsubjpath(aap,subj),['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_design.jpg']));
            close all;
        end
        fdiag = dir(fullfile(aas_getsubjpath(aap,subj),'diagnostic_*.jpg'));
        for d = 1:numel(fdiag)
            aap = aas_report_add(aap,subj,'<table><tr><td>');
            aap=aas_report_addimage(aap,subj,fullfile(aas_getsubjpath(aap,subj),fdiag(d).name));
            aap = aas_report_add(aap,subj,'</td></tr></table>');
        end

    case 'doit'
        % Get subject directory
        cwd=pwd;

        % We can now have missing sessions per subject, so we're going to use only
        % the sessions that are common to this subject and selected_sessions
        [numSess, sessInds] = aas_getN_bydomain(aap, 'session', subj);
        subjSessionI = intersect(sessInds, aap.acq_details.selected_sessions);
        
        % Prepare basic SPM model...
        [SPM, anadir, files, allfiles, model, modelC] = aas_firstlevel_model_prepare(aap, subj);
        
        % Get all the nuisance regressors...
        [movementRegs, compartmentRegs, physiologicalRegs, spikeRegs, GLMDNregs] = ...
            aas_firstlevel_model_nuisance(aap, subj, files);

        %% Set up CORE model
        cols_nuisance=[];
        cols_interest=[];
        currcol=1;

        sessnuminspm=0;

        for sess = 1:numSess
            sessnuminspm=sessnuminspm+1;

            % Settings
            SPM.nscan(sessnuminspm) = size(files{sess},1);
            SPM.xX.K(sessnuminspm).HParam = aap.tasklist.currenttask.settings.highpassfilter;

            % Set up model
            [SPM, cols_interest, cols_nuisance, currcol] = ...
                aas_firstlevel_model_define(aap, sess, sessnuminspm, SPM, model, modelC, ...
                                                             cols_interest, cols_nuisance, currcol, ...
                                                             movementRegs, compartmentRegs, physiologicalRegs, spikeRegs, GLMDNregs);
        end

        cd (anadir)

        %%%%%%%%%%%%%%%%%%%
        %% DESIGN MATRIX %%
        %%%%%%%%%%%%%%%%%%%
        SPM.xY.P = allfiles;
        SPMdes = spm_fmri_spm_ui(SPM);

        SPMdes.xX.X = double(SPMdes.xX.X);
        
        % DIAGNOSTIC
        mriname = aas_prepare_diagnostic(aap, subj);
        try
            saveas(1, fullfile(aap.acq_details.root, 'diagnostics', ...
                                                [mfilename '__' mriname '.fig']));
        catch
        end

        % now check real covariates and nuisance variables are
        % specified correctly
        SPMdes.xX.iG = cols_nuisance;
        SPMdes.xX.iC = cols_interest;

        % Turn off masking if requested
        if ~aap.tasklist.currenttask.settings.firstlevelmasking
            SPMdes.xM.I=0;
            SPMdes.xM.TH=-inf(size(SPMdes.xM.TH));
        end

        %%%%%%%%%%%%%%%%%%%
        %% ESTIMATE MODEL%%
        %%%%%%%%%%%%%%%%%%%
        spm_unlink(fullfile('.', 'mask.img')); % avoid overwrite dialog
        SPMest = spm_spm(SPMdes);

        %% Describe outputs
        cd (cwd);

        % Describe outputs
        %  firstlevel_spm
        aap=aas_desc_outputs(aap,subj,'firstlevel_spm',fullfile(anadir,'SPM.mat'));

        %  firstlevel_betas (includes related statistical files)
        allbetas=dir(fullfile(anadir,'beta_*'));
        allbetas=vertcat(allbetas,...
            dir(fullfile(anadir,'ResMS.*')),...
            dir(fullfile(anadir,'RPV.*')));
        if aap.tasklist.currenttask.settings.firstlevelmasking
            allbetas=vertcat(allbetas,...
                dir(fullfile(anadir,'mask.*')));
        end
        betafns=strcat(repmat([anadir filesep],[numel(allbetas) 1]),char({allbetas.name}));
        aap=aas_desc_outputs(aap,subj,'firstlevel_betas',betafns);

        %% DIAGNOSTICS...
%         h = firstlevelmodelStats(anadir, [], fullfile(anadir, 'mask.img'));
%         saveas(h.regs, fullfile(aap.acq_details.root, 'diagnostics', [mfilename '__' mriname '_regs.eps']), 'psc2');
%         saveas(h.betas, fullfile(aap.acq_details.root, 'diagnostics', [mfilename '__' mriname '_betas.eps']), 'psc2');
%         
%         close(h.regs)
%         close(h.betas)
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end