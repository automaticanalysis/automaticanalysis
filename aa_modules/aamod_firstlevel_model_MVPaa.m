% AA module - first level statistics
% This module is based on the work by Mumford, Turner, Ashby & Poldrack,
% 2012, Neuroimage 59 (2012)
% Essentially, this makes a separate model for each regressor, allowing to
% model events such that the betas extracted are less correlated with one
% another, theoretically improving power in the MVPA analysis
% Meant to be used for rapid event-related designs
% **********************************************************************
% You should no longer need to change this module - you may just
% modify the .xml or model in your user script
% **********************************************************************
% Based on original by FIL London and Adam Hampshire MRC CBU Cambridge Feb 2006
% Modified for aa by Rhodri Cusack MRC CBU Mar 2006-Aug 2007
% Thanks to Rik Henson for various suggestions

function [aap,resp]=aamod_firstlevel_model_MVPaa(aap,task,subj)

resp='';

switch task
    case 'report'
        
    case 'doit'
        %get subject directory
        cwd=pwd;
        % This does not work [AVG]
        %[junk, subjname]=fileparts(subj_dir);
        % Try this instead!
        subjname = aap.acq_details.subjects(subj).mriname;
        
        %% Movement regressors (extended!) [AVG]
        [moves, mnames] = aas_movPars(aap,subj, aap.tasklist.currenttask.settings.moveMat);
        
        %% Compartment regressors [AVG]
        compRegNames = {'GM', 'WM', 'CSF', 'OOH'};
        compTC = [];
        Cregs = cell(1,length(aap.acq_details.sessions));
        for sess=aap.acq_details.selected_sessions
            % If we don't have compartment Signals, this should give up...
            try
                load(aas_getfiles_bystream(aap,subj,sess,'compSignal'));
                Cregs{sess} = compTC;
            catch
            end
        end
        
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
        % TR is manually specified (not recommended as source of error)
        if isfield(aap.tasklist.currenttask.settings,'TR') && ...
                ~isempty(aap.tasklist.currenttask.settings.TR)
            SPM.xY.RT =aap.tasklist.currenttask.settings.TR;
        else
            % Get TR from DICOM header checking they're the same for all sessions
            for sess=aap.acq_details.selected_sessions
                DICOMHEADERS=load(aas_getfiles_bystream(aap,subj,sess,'epi_dicom_header'));
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
        % can have module specific value, but kept for backwards compatability
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
        
        %% Set up CORE model
        coreSPM = SPM;
        cols_nuisance=[];
        cols_interest=[];
        sessnuminspm=1;
        currcol=1;
        
        % [AVG] - let's save the model...
        model = cell(size(aap.acq_details.sessions));
        files = cell(size(aap.acq_details.sessions));
        
        for sess = aap.acq_details.selected_sessions
            % Get model{sess} data from aap
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
            
            %% Should now have just one model{sess} spec
            modelnum=find(sessmatches & subjmatches);
            if (length(modelnum)>1)
                aas_log(aap,true,sprintf('Error while getting model details as more than one specification for subject %s session %s',subjname,aap.acq_details.sessions(sess).name));
            end
            if (isempty(modelnum))
                aas_log(aap,true,'Cannot find model specification. Check either user script or aamod_firstlevel_model{sess}.xml');
            end
            
            model{sess}=aap.tasklist.currenttask.settings.model(modelnum);
                        
            files{sess} = aas_getimages_bystream(aap,subj,sess,'epi');
            
            SPM.nscan(sessnuminspm) = size(files{sess},1);
            
            allfiles = strvcat(allfiles,files{sess});
            
            for c = 1:length(model{sess}.event);
                if (isempty(model{sess}.event(c).parametric))
                    parametric=struct('name','none');
                else
                    parametric=model{sess}.event(c).parametric;
                end
                SPM.Sess(sessnuminspm).U(c) = struct(...
                    'ons',model{sess}.event(c).ons,...
                    'dur',model{sess}.event(c).dur,...
                    'name',{{model{sess}.event(c).name}},...
                    'P',parametric);
                cols_interest=[cols_interest currcol];
                currcol=currcol+1;
            end
            
            SPM.xX.K(sessnuminspm).HParam = aap.tasklist.currenttask.settings.highpassfilter;
            
            %% Movement and other nuisance regressors: compartments [AVG]
            if aap.tasklist.currenttask.settings.includemovementpars==1
                SPM.Sess(sessnuminspm).C.C    = [moves{sess} Cregs{sess}(:, aap.tasklist.currenttask.settings.compRegs)];     % [n x c double] covariates
                SPM.Sess(sessnuminspm).C.name = [mnames compRegNames{aap.tasklist.currenttask.settings.compRegs}]; % [1 x c cell]   names
                cols_nuisance=[cols_nuisance (currcol:(currcol+length(mnames)+length(aap.tasklist.currenttask.settings.compRegs)-1))];
                currcol=currcol+length(mnames)+length(aap.tasklist.currenttask.settings.compRegs);
            else
                SPM.Sess(sessnuminspm).C.C = [];
                SPM.Sess(sessnuminspm).C.name = {};
            end
            sessnuminspm=sessnuminspm+1;
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
        
        %% REDO MODEL WITH MVPA bonus...
        
        % Make temporary directory inside folder
        Tanadir = (fullfile(anadir, 'temp'));
        
        % Find out how large the model{sess} should be (per session)
        sessRegs = 1:size(model{aap.acq_details.selected_sessions(1)}.event,2);
        
        for n = sessRegs
            try rmdir(Tanadir); catch; end
            mkdir(Tanadir)
        
            cols_nuisance=[];
            cols_interest=[];
            sessnuminspm=1;
            currcol=1;
            rSPM = coreSPM;
        
            for sess = aap.acq_details.selected_sessions
                % * Get noise regressor numbers
                noiseRegs = sessRegs;
                noiseRegs(n) = [];
                
                % No need to check model{sess} again, did already before...
                Tmodel = model{sess};
                Tmodel.event = [];
                % Event names, rather simple...
                Tmodel.event(1).name = 'Reg';
                Tmodel.event(2).name = 'Noise';
                % I don't see a way to do parameteric stuff here
                Tmodel.event(1).parametric = [];
                Tmodel.event(2).parametric = [];
                % Easy to set up the Regressor event
                Tmodel.event(1).ons = model{sess}.event(n).ons;
                Tmodel.event(1).dur = model{sess}.event(n).dur;
                % Trickier for the Noise event
                Tons = [];
                Tdur = [];
                for r = sessRegs
                    % Make sure we don't include the modelled regressor in
                    % the noise regressor
                    if r ~= n
                        Tons = [Tons; model{sess}.event(r).ons];
                        Tdur = [Tdur; model{sess}.event(r).dur];
                    end
                end
                
                % Sort the onsets, and apply same reordering to durations
                [Tons, ind]=sort(Tons);
                if (length(Tdur)>1)
                    Tdur=Tdur(ind);
                end
                Tmodel.event(2).ons = Tons;
                Tmodel.event(2).dur = Tdur;
                
                rSPM.nscan(sessnuminspm) = size(files{sess},1);
                
                for c = 1:length(Tmodel.event);
                    if (isempty(model{sess}.event(c).parametric))
                        parametric=struct('name','none');
                    else
                        parametric=Tmodel.event(c).parametric;
                    end
                    rSPM.Sess(sessnuminspm).U(c) = struct(...
                        'ons',Tmodel.event(c).ons,...
                        'dur',Tmodel.event(c).dur,...
                        'name',{{Tmodel.event(c).name}},...
                        'P',parametric);
                    cols_interest=[cols_interest currcol];
                    currcol=currcol+1;
                end
                
                rSPM.xX.K(sessnuminspm).HParam = aap.tasklist.currenttask.settings.highpassfilter;
                
                %% Movement and other nuisance regressors: compartments [AVG]
                if aap.tasklist.currenttask.settings.includemovementpars==1
                    rSPM.Sess(sessnuminspm).C.C    = [moves{sess} Cregs{sess}(:, aap.tasklist.currenttask.settings.compRegs)];     % [n x c double] covariates
                    rSPM.Sess(sessnuminspm).C.name = [mnames compRegNames{aap.tasklist.currenttask.settings.compRegs}]; % [1 x c cell]   names
                    cols_nuisance=[cols_nuisance (currcol:(currcol+length(mnames)+length(aap.tasklist.currenttask.settings.compRegs)-1))];
                    currcol=currcol+length(mnames)+length(aap.tasklist.currenttask.settings.compRegs);
                else
                    rSPM.Sess(sessnuminspm).C.C = [];
                    rSPM.Sess(sessnuminspm).C.name = {};
                end
                sessnuminspm=sessnuminspm+1;
            end
            cd(Tanadir)

            % DEBUG
            %{
            if n == 1
                subplot(3,1,1); plot(SPMdes.xX.X(1:100,1:sessRegs(end))); title('Normal model')
                subplot(3,1,2); plot(rSPMdes.xX.X(1:100,1:2)); title('New model')
                subplot(3,1,3); plot([SPMdes.xX.X(1:100,1) sum(SPMdes.xX.X(1:100,2:sessRegs(end)),2)]); title('Old model, summed')
            end
            %}
            
            rSPM.xY.P = allfiles;
            rSPMdes = spm_fmri_spm_ui(rSPM);
            
            % now check real covariates and nuisance variables are
            % specified correctly
            rSPMdes.xX.iG=cols_nuisance;
            rSPMdes.xX.iC=cols_interest;
            
            spm_unlink(fullfile('.', 'mask.img')); % avoid overwrite dialog
            rSPMest = spm_spm(rSPMdes);
            
            % After this, move the betas to correct location
            % Find the regressors we wish to move...
            % odd ones in rSPMest, since there's only 2 regressors...
            nOrig = rSPMest.xX.iC(1: ... 
                2 ... 
                :length(rSPMest.xX.iC));
            % determined by the number of regressors per session in SPMest
            nDest = SPMest.xX.iC((0: ...
                length(sessRegs): ...
                length(sessRegs)*length(aap.acq_details.selected_sessions)-1)+n);
            
            % Now move the actual files
            for f = 1:length(nOrig)
                unix(['mv ' fullfile(Tanadir, sprintf('beta_%04d.img', nOrig(f))) ...
                    ' ' fullfile(anadir, sprintf('beta_%04d.img', nDest(f)))]);
                unix(['mv ' fullfile(Tanadir, sprintf('beta_%04d.hdr', nOrig(f))) ...
                    ' ' fullfile(anadir, sprintf('beta_%04d.hdr', nDest(f)))]);
            end
        end
        try rmdir(Tanadir); catch; end
        
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