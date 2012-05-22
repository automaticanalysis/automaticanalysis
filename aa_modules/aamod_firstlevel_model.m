% AA module - first level statistics
% **********************************************************************
% You should no longer need to change this module - you may just
% modify the .xml or model in your user script
% **********************************************************************
% Based on original by FIL London and Adam Hampshire MRC CBU Cambridge Feb 2006
% Modified for aa by Rhodri Cusack MRC CBU Mar 2006-Aug 2007
% Thanks to Rik Henson for various suggestions

function [aap,resp]=aamod_firstlevel_model(aap,task,subj)

resp='';

switch task
    case 'report'
        
    case 'doit'
        %get subject directory
        cwd=pwd;
        % This does not work [AVG]
        %[~, subjname]=fileparts(subj_dir);
        % Try this instead!
        subjname = aap.acq_details.subjects(subj).mriname;
        
        %% Movement regressors (extended!) [AVG]
        [moves, mnames] = aas_movPars(aap,subj, aap.tasklist.currenttask.settings.moveMat);
        
        %% Compartment regressors [AVG]
        compRegNames = {'GM', 'WM', 'CSF', 'OOH'};
        compTC = [];
        Cregs = cell(1,length(aap.acq_details.sessions));
        for sess = aap.acq_details.selected_sessions
            % If we have compartment Signals load them...
            if (exist(aas_getinputstreamfilename(aap,subj,sess),'file'))
                load(aas_getfiles_bystream(aap,subj,sess,'compSignal'));
            end
        end
        
        % Covariates that may or may not exist [AVG]
        covars = cell(1,length(aap.acq_details.sessions));
        
        clear SPM
        
        %% Set up basis functions
        if (isfield(aap.tasklist.currenttask.settings,'xBF'))
            SPM.xBF=aap.tasklist.currenttask.settings.xBF;
        else
            SPM.xBF.T          = 16; % number of time bins per scan
            SPM.xBF.UNITS      = 'scans';           % OPTIONS: 'scans'|'secs' for onsets
            SPM.xBF.Volterra   = 1;                 % OPTIONS: 1|2 = order of convolution
            SPM.xBF.name       = 'hrf';
            SPM.xBF.length     = 32;                % length in seconds
            SPM.xBF.order      = 1;                 % order of basis set
        end
        
        firstsess=aap.acq_details.selected_sessions(1);
        
        %% retrieve TR from DICOM header
        % if TR is manually specified (not recommended as source of error)
        if isfield(aap.tasklist.currenttask.settings,'TR') && ...
                ~isempty(aap.tasklist.currenttask.settings.TR)
            SPM.xY.RT =aap.tasklist.currenttask.settings.TR;
        else
            % Get TR from DICOM header checking they're the same for all sessions
            for sess=aap.acq_details.selected_sessions
                DICOMHEADERS=load(aas_getfiles_bystream(aap,subj,sess,'epi_header'));
                try
                    TR=DICOMHEADERS.DICOMHEADERS{1}.volumeTR;
                catch
                    % [AVG] This is for backwards compatibility!
                    TR=DICOMHEADERS.DICOMHEADERS{1}.RepetitionTime/1000;
                end
                if (sess==firstsess)
                    SPM.xY.RT = TR;
                else
                    if (SPM.xY.RT~=TR)
                        aas_log(aap,true,sprintf('Session %d has different TR from earlier sessions, they can''t be in the same model.',sess));
                    end
                end
            end
        end
        
        %% Get slice order from sliceorder stream if it exists, check same
        % for all sessions
        usesliceorder=aas_stream_has_contents(aap,'sliceorder');
        if (usesliceorder)
            for sess=aap.acq_details.selected_sessions
                sliceorderstruct=load(aas_getfiles_bystream(aap,subj,sess,'sliceorder'));
                if (sess==firstsess)
                    sliceorder=sliceorderstruct.sliceorder;
                    refslice=sliceorderstruct.refslice;
                else
                    if (any(sliceorderstruct.sliceorder~=sliceorder))
                        aas_log(aap,true,sprintf('Session %d has different slice order from earlier sessions, they can''t be in the same model.',sess));
                    end
                end
            end
        end
        
        SPM.xGX.iGXcalc = 'None';
        SPM.xVi.form = 'AR(1)';
        
        %% Adjust time bin T0 according to reference slice & slice order
        %  implements email to CBU from Rik Henson 27/06/07
        %  assumes timings are relative to beginning of scans
        
        if (usesliceorder)
            refwhen=(find(sliceorder==refslice)-1)/(length(sliceorder)-1);
        else
            aas_log(aap,false,'No stream sliceorder found, defaulting timing to SPM.xBF.T0=0 in model');
            refwhen=0;
        end
        SPM.xBF.T0 = round(SPM.xBF.T*refwhen);
        
        subdata = aas_getsubjpath(aap,subj);
        
        %% Deal with extraparameters. Not needed any more, as
        % aap.directory_conventions.stats_singlesubj
        % can have module specific value, but kept for backwards
        % compatability
        if (isfield(aap.tasklist.currenttask.extraparameters,'stats_suffix'))
            stats_suffix=aap.tasklist.currenttask.extraparameters.stats_suffix;
        else
            stats_suffix=[];
        end
        
        anadir = fullfile(subdata,[aap.directory_conventions.stats_singlesubj stats_suffix]);
        if ~exist(anadir,'dir')
            mkdir(subdata,[aap.directory_conventions.stats_singlesubj stats_suffix]);
        end
        
        allfiles='';
        
        %% Set up model
        cols_nuisance=[];
        cols_interest=[];
        sessnuminspm=0;
        currcol=1;
        for sess = aap.acq_details.selected_sessions
            sessnuminspm=sessnuminspm+1;
            
            %% Get model data from aap
            subjmatches=strcmp(subjname,{aap.tasklist.currenttask.settings.model.subject});
            sessmatches=strcmp(aap.acq_details.sessions(sess).name,{aap.tasklist.currenttask.settings.model.session});
            % If no exact spec found, try session wildcard, then subject
            % wildcard, then wildcard for both
            if (~any(sessmatches & subjmatches))
                sesswild=strcmp('*',{aap.tasklist.currenttask.settings.model.session});
                if (any(sesswild & subjmatches))
                    sessmatches=sesswild;
                else
                    subjwild=strcmp('*',{aap.tasklist.currenttask.settings.model.subject});
                    if (any(sessmatches & subjwild))
                        subjmatches=subjwild;
                    else
                        subjmatches=subjwild;
                        sessmatches=sesswild;
                    end
                end
            end
            
            % Should now have just one model spec
            modelnum=find(sessmatches & subjmatches);
            
            %% Get modelC (covariate) data from aap
            subjmatches=strcmp(subjname,{aap.tasklist.currenttask.settings.modelC.subject});
            sessmatches=strcmp(aap.acq_details.sessions(sess).name,{aap.tasklist.currenttask.settings.modelC.session});
            % If no exact spec found, try session wildcard, then subject
            % wildcard, then wildcard for both
            if (~any(sessmatches & subjmatches))
                sesswild=strcmp('*',{aap.tasklist.currenttask.settings.modelC.session});
                if (any(sesswild & subjmatches))
                    sessmatches=sesswild;
                else
                    subjwild=strcmp('*',{aap.tasklist.currenttask.settings.modelC.subject});
                    if (any(sessmatches & subjwild))
                        subjmatches=subjwild;
                    else
                        subjmatches=subjwild;
                        sessmatches=sesswild;
                    end
                end
            end
            
            % Should now have just one modelC spec
            modelCnum=find(sessmatches & subjmatches);
            
            %% Check that we have at least one model of interest (normal or covariate)
            if (length(modelnum)>1) || (length(modelCnum)>1)
                aas_log(aap,true,sprintf('Error while getting model details as more than one specification for subject %s session %s',subjname,aap.acq_details.sessions(sess).name));
            end
            if (isempty(modelnum)) && (isempty(modelCnum))
                aas_log(aap,true,'Cannot find model specification. Check either user script or aamod_firstlevel_model.xml');
            end
            
            %% Now do the model...
            
            if ~isempty(modelnum)
                model=aap.tasklist.currenttask.settings.model(modelnum);
                for c = 1:length(model.event);
                    if (isempty(model.event(c).parametric))
                        parametric=struct('name','none');
                    else
                        parametric=model.event(c).parametric;
                    end
                    SPM.Sess(sessnuminspm).U(c) = struct(...
                        'ons',model.event(c).ons,...
                        'dur',model.event(c).dur,...
                        'name',{{model.event(c).name}},...
                        'P',parametric);
                    cols_interest=[cols_interest currcol];
                    currcol=currcol+1;
                end
            else
                SPM.Sess(sessnuminspm).U = [];
                SPM.Sess(sessnuminspm).row = [];
            end
            
            %% ADD ANY EXISTING COVARIATES HERE...
            
            SPM.Sess(sessnuminspm).C.C = [];
            SPM.Sess(sessnuminspm).C.name = {};
            
            if ~isempty(modelCnum)
                
                modelC = aap.tasklist.currenttask.settings.modelC(modelCnum);
                
                %% Set up the convolution vector...
                % xBF.dt      - time bin length {seconds}
                % xBF.name    - description of basis functions specified
                % xBF.length  - window length (seconds)
                % xBF.order   - order
                
                xBF = [];
                xBF.dt = SPM.xY.RT;
                xBF.name = SPM.xBF.name;
                xBF.length = SPM.xBF.length;
                xBF.order = SPM.xBF.order;
                xBF = spm_get_bf(xBF);
                
                for c = 1:length(modelC.covariate);
                    covVect = modelC.covariate(c).vector;
                    
                    % Do we convolve with HRF?
                    if modelC.covariate(c).HRF > 0
                        U =[];
                        U.u = covVect;
                        U.name = {modelC.covariate(c).name};
                        covVect = spm_Volterra(U, xBF.bf);
                    end
                    
                    SPM.Sess(sessnuminspm).C.C    = [SPM.Sess(sessnuminspm).C.C ...
                        covVect];     % covariate
                    SPM.Sess(sessnuminspm).C.name = [SPM.Sess(sessnuminspm).C.name ...
                        modelC.covariate(c).name];
                    
                    % Is the covariate of interest or nuisance
                    if modelC.covariate(c).interest > 0
                        cols_interest=[cols_interest currcol];
                    else
                        cols_nuisance=[cols_nuisance currcol];
                    end
                    currcol=currcol + 1;
                end
            end
            
            %% Movement and other nuisance regressors: compartments [AVG]
            if aap.tasklist.currenttask.settings.includemovementpars==1
                SPM.Sess(sessnuminspm).C.C    = [SPM.Sess(sessnuminspm).C.C ...
                    moves{sess} ...
                    Cregs{sess}(:, aap.tasklist.currenttask.settings.compRegs)];     % [n x c double] covariates
                SPM.Sess(sessnuminspm).C.name = [SPM.Sess(sessnuminspm).C.name ...
                    mnames ...
                    compRegNames{aap.tasklist.currenttask.settings.compRegs}]; % [1 x c cell]   names
                cols_nuisance=[cols_nuisance [currcol:(currcol + ...
                    length(mnames) ...
                    + length(aap.tasklist.currenttask.settings.compRegs) - 1)]];
                currcol=currcol ...
                    + length(mnames) ...
                    + length(aap.tasklist.currenttask.settings.compRegs);
            end
            
            %% SETTINGS & GET FILES
            
            files = aas_getfiles_bystream(aap,subj,sess,'epi');
            allfiles = strvcat(allfiles,files);
            
            SPM.xX.K(sessnuminspm).HParam = aap.tasklist.currenttask.settings.highpassfilter;
            
            SPM.nscan(sessnuminspm) = size(files,1);
        end
        
        cd (anadir)
        
        SPM.xY.P = allfiles;
        SPMdes = spm_fmri_spm_ui(SPM);
        
        % now check real covariates and nuisance variables are
        % specified correctly
        SPMdes.xX.iG=cols_nuisance;
        SPMdes.xX.iC=cols_interest;
        
        spm_unlink(fullfile('.', 'mask.img')); % avoid overwrite dialog
        SPMest = spm_spm(SPMdes);
        
        %% Describe outputs
        
        cd (cwd);
        
        % Describe outputs
        %  firstlevel_spm
        aap=aas_desc_outputs(aap,subj,'firstlevel_spm',fullfile(anadir,'SPM.mat'));
        
        %  firstlevel_betas (includes related statistical files)
        allbetas=dir(fullfile(anadir,'beta_*'));
        betafns=[];
        for betaind=1:length(allbetas);
            betafns=strvcat(betafns,fullfile(anadir,allbetas(betaind).name));
        end
        otherfiles={'mask.hdr','mask.img','ResMS.hdr','ResMS.img','RPV.hdr','RPV.img'};
        for otherind=1:length(otherfiles)
            betafns=strvcat(betafns,fullfile(anadir,otherfiles{otherind}));
        end
        aap=aas_desc_outputs(aap,subj,'firstlevel_betas',betafns);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end