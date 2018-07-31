% AA module - first level statistics 1: configure model (add events and
% covariates)
%
% Split version of aamod_firstlevel_model - useful for inserting custom model
% options in between stages. Replaces aamod_firstlevel_model as a whole. Stages:
%
% aamod_firstlevel_model_1_config - add events, set up basic SPM struct
% aamod_firstlevel_model_2_convolve - generate convolved design matrix
% (spm_fmri_spm_ui)
% aamod_firstlevel_model_3_estimate - fit the model (spm_spm)

function [aap,resp]=aamod_firstlevel_model_1_config(aap,task,subj)

resp='';

switch task
    case 'report' 
    case 'doit'
        % Get subject directory
        cwd=pwd;
        
 		% We can now have missing sessions per subject, so we're going to use only
        % the sessions that are common to this subject and selected_sessions
        [numSess, sessInds] = aas_getN_bydomain(aap, 'session', subj);
        subjSessionI = intersect(sessInds, aap.acq_details.selected_sessions);
        numSess = numel(subjSessionI);

        % Add PPI if exist
        modname = aap.tasklist.currenttask.name;
        modnameind = regexp(modname, '_\d{5,5}$');
        modindex = str2num(modname(modnameind+1:end));
        for sess = subjSessionI
            if aas_stream_has_contents(aap,subj,sess,'ppi')
                load(aas_getfiles_bystream(aap,subj,sess,'ppi'));
                [phys, psych] = strtok(PPI.name,'x('); psych = psych(3:end-1);
                aap = aas_addcovariate(aap,modname,...
                    basename(aas_getsubjpath(aap,subj)),aap.acq_details.sessions(sess).name,...
                    'PPI',PPI.ppi,0,1);
                aap = aas_addcovariate(aap,modname,...
                    basename(aas_getsubjpath(aap,subj)),aap.acq_details.sessions(sess).name,...
                    ['Psych_' psych],PPI.P,0,1);
                aap = aas_addcovariate(aap,modname,...
                    basename(aas_getsubjpath(aap,subj)),aap.acq_details.sessions(sess).name,...
                    ['Phys_' phys],PPI.Y,0,1);
            end
        end
        % update current sesstings
        aap.tasklist.currenttask.settings.modelC = aap.tasksettings.(modname(1:modnameind-1))(modindex).modelC;
        
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
        % these fields are only set after convolving, but we specify them now so
        % we can access them in aamod_firstlevel_model_2_convolve.
        SPM.xX.iG = cols_nuisance;
        SPM.xX.iC = cols_interest;
        SPM.xY.P = allfiles;

        % Describe outputs
        %  firstlevel_spm
        outfile = fullfile(anadir,'SPM.mat');
        save(outfile,'SPM');
        aap=aas_desc_outputs(aap,subj,'firstlevel_spm',outfile);
        cd (cwd);

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
