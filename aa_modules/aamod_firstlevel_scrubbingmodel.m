function [aap,resp]=aamod_firstlevel_scrubbingmodel(aap, task, subjInd)
%
% This module creates a first level model intended to denoise fMRI
% timeseries.  Columns of the GLM can include: Bandpass filter regressors,
% CSF and WM signals, motion parameters, their lag-3 2nd order volteraa 
% expansion, and "spike" regressors".  If you don't want to include any of
% those things, set the apporporiate option in the corresponding .xml file,
% and comment out any streams you don't need.
%
% It should also be possible to include task-related regressors if you wish
% to regress those out the data, too.
%
% Once the model is created, aamod_firstlevel_modelestimate_saveresids
% should be rrun to actually estimate the model and save the residual time
% series.
%
% This module is based on one originally created by Jonathan Peelle, and
% includes snippets of code supplied by Rik Henson.  Thank you!
%
% -- conorwild -- 02/12/2013 

resp='';

switch task
    case 'report'
        
    case 'doit'
        
        settings = aap.tasklist.currenttask.settings;
        
        startingDir = pwd();
        
        % Get options
        highpassFilter = settings.highpassfilter;
        if isempty(highpassFilter), highpassFilter = Inf; end
        
        TR = settings.TR;
        globalScaling = settings.globalscaling;
        autocorrelation = settings.autocorrelation;
        includeMovement = settings.includemovementpars;   % include movement parameters?
        volterraMovement = settings.volterramovementpars; % include Volterra expansion of movement parameters? (passed to aas_movPars)
        bandPass = settings.bandpass;                     % for resting state functional connectivity
        svdThresh = settings.svdthresh;                   % if < 1, do dimension reduction on confounds
        explicitMask = settings.explicitmask;
        maskThreshold = settings.maskthreshold;
        
        subjName = aap.acq_details.subjects(subjInd).mriname;
        subjPath = aas_getsubjpath(aap, subjInd);
        
        % Deal with extraparameters. Not needed any more, as
        % aap.directory_conventions.stats_singlesubj
        % can have module specific value, but kept for backwards
        % compatability
        if (isfield(aap.tasklist.currenttask.extraparameters,'stats_suffix'))
            stats_suffix=aap.tasklist.currenttask.extraparameters.stats_suffix;
        else
            stats_suffix=[];
        end
        
        analysisDir = fullfile(subjPath,[aap.directory_conventions.stats_singlesubj stats_suffix]);
        if ~exist(analysisDir,'dir')
            mkdir(analysisDir);
        end
        cd(analysisDir);
        
        SPM = [];
        SPM.swd = analysisDir;
        
        global defaults
        defaults.mask.thresh = maskThreshold;

        % Get basis functions from task settings
        SPM.xBF = settings.xBF;
        
        firstsess = aap.acq_details.selected_sessions(1);
        
        % if TR is manually specified, use that, otherwise try to get from DICOMs.
        if isfield(settings,'TR') && ~isempty(TR)
            SPM.xY.RT = TR;
        else
            try
                
                % Get TR from DICOM header checking they're the same for all sessions
                
                for sess = aap.acq_details.selected_sessions
                    DICOMHEADERS=load(aas_getfiles_bystream(aap, subjInd, sess, 'epi_dicom_header'));
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
                            aas_log(aap, true, sprintf('Session %d has different TR from earlier sessions, they can''t be in the same model.', sess));
                        end
                    end
                end
            catch
                aas_log(aap, true, 'No epi_dicom_header information available - please set TR yourself');
            end
        end
        
        
        SPM.xGX.iGXcalc = globalScaling; % 'None';
        SPM.xVi.form = autocorrelation; % 'AR(1)';    
        
        % Loop through sessions and specify model
        nuisanceCols = [];
        interestCols = [];
        spmSession = 1; % which session in SPM
        currentCol = 1;
        allFiles = [];
        
        for thisSess = aap.acq_details.selected_sessions       
            
            % Any session-specific settings?
            SPM.xX.K(spmSession).HParam = highpassFilter;
            
            % Get files (do this first as some other options need to know
            % how many files)
            thisSessFiles = aas_getfiles_bystream(aap, subjInd, thisSess, 'epi');
            nScans = size(thisSessFiles, 1);      
            
            % Get model data from aap
            subjmatches=strcmp(subjName,{settings.model.subject});
            sessmatches=strcmp(aap.acq_details.sessions(thisSess).name,{settings.model.session});
            % If no exact spec found, try session wildcard, then subject
            % wildcard, then wildcard for both
            if (~any(sessmatches & subjmatches))
                sesswild=strcmp('*',{settings.model.session});
                if (any(sesswild & subjmatches))
                    sessmatches=sesswild;
                else
                    subjwild=strcmp('*',{settings.model.subject});
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
            
            % Get modelC (covariate) data from aap
            subjmatches=strcmp(subjName,{settings.modelC.subject});
            sessmatches=strcmp(aap.acq_details.sessions(thisSess).name,{settings.modelC.session});
            
            % If no exact spec found, try session wildcard, then subject
            % wildcard, then wildcard for both
            if (~any(sessmatches & subjmatches))
                sesswild=strcmp('*',{settings.modelC.session});
                if (any(sesswild & subjmatches))
                    sessmatches=sesswild;
                else
                    subjwild=strcmp('*',{settings.modelC.subject});
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
            
            % Check that we have at least one model of interest (normal or covariate)
            if (length(modelnum)>1) || (length(modelCnum)>1)
                aas_log(aap,true,sprintf('Error while getting model details as more than one specification for subject %s session %s',subjName,aap.acq_details.sessions(sess).name));
            end
            
            if ~isempty(modelnum)
                model=settings.model(modelnum);
                for c = 1:length(model.event);
                    if (isempty(model.event(c).parametric))
                        parametric=struct('name','none');
                    else
                        parametric=model.event(c).parametric;
                    end
                    SPM.Sess(spmSession).U(c) = struct(...
                        'ons',model.event(c).ons,...
                        'dur',model.event(c).dur,...
                        'name',{{model.event(c).name}},...
                        'P',parametric);
                    cols_interest=[cols_interest currcol];
                    currcol=currcol+1;
                end
            else
                SPM.Sess(spmSession).U = [];
                SPM.Sess(spmSession).row = [];
            end
            
            % Covariates
            SPM.Sess(spmSession).C.C = [];
            SPM.Sess(spmSession).C.name = {};
            
            if ~isempty(modelCnum)
                
                modelC = settings.modelC(modelCnum);
                
                % Set up the convolution vector...
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
                        U.u = covVect(:);
                        U.name = {modelC.covariate(c).name};
                        covVect = spm_Volterra(U, xBF.bf);
                    end
                    
                    SPM.Sess(spmSession).C.C    = [SPM.Sess(spmSession).C.C ...
                        covVect];     % covariate
                    SPM.Sess(spmSession).C.name = [SPM.Sess(spmSession).C.name ...
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

            % Confounds/covariates of no interest/whatever you want to call
            % them. These are collected together in case one wants to do
            % dimensionality reduction.
            
            C = []; % confounds
            Cnames = []; % names
          
            % Compartment signals (GM, WM, CSF)
            compFiles = aas_getfiles_bystream(aap, subjInd, spmSession, 'compSignal');
            compSignals = load(compFiles);
            
            if size(compSignals.compTC, 1) ~= nScans, aas_log(aap, 'compSig is wrong length!', 1); end
            
            % Include a CSF signal?
            if settings.includeCSF
                fprintf('Adding CSF signal as covariate.\n');
                C = [C compSignals.compTC(:,3)];
                Cnames = [Cnames 'CSF'];
            end
            
            % Include a WM signal?
            if settings.includeCSF
                fprintf('Adding WM signal as covariate.\n');
                C = [C compSignals.compTC(:,2)];
                Cnames = [Cnames 'WM'];
            end
            
            % Add movement parameters? (this includes Volterra expansion, if requested - all gotten before looping through sessions)
            if includeMovement
                
                M = spm_load(aas_getfiles_bystream(aap, subjInd, spmSession, 'realignment_parameter'));
                
                if volterraMovement
                    U=[];
                    for c = 1 : 6
                        U(c).u = M(:,c); 
                        U(c).name{1}='c'; 
                    end
                    M = spm_Volterra(U, [1 0 0; 1 -1 0; 0 1 -1]', 2);
                end
                
                C = [C M];
                Cnames = [Cnames arrayfun(@(x) sprintf('Mov%d', x), [1:size(M, 2)], 'UniformOutput', false)];
                aas_log(aap, false, sprintf('%d movement regressors added for this session.', size(M,2)));
            end
            
            % Bandpass filtering (using DCT)
            if isempty(bandPass)
                numColK = 0;
            else
                
                if ischar(bandPass)
                    bandPass = str2num(bandPass); % in case specified as a string
                end
                
                highPassCut = 1/bandPass(1);
                lowPassCut = 1/bandPass(2); % They were specified in Hz!
   
                K = spm_dctmtx(nScans, nScans);
                nHP = fix(2*(nScans*TR)/highPassCut + 1);
                nLP = fix(2*(nScans*TR)/lowPassCut + 1);
                
                K = K(:, [1:nHP nLP:nScans]); % includes constant
                numColK  = size(K,2);
                
                C = [C K];
                
                Cnames = [Cnames arrayfun(@(x) sprintf('BP%d', x), [1:size(K, 2)], 'UniformOutput', false)];
                
                aas_log(aap, false, sprintf('%d temporal filtering regressors added for this session.', size(K,2)));
            end
            
            if settings.includeSpikes
                
                % Get the spikes and moves
                load(aas_getfiles_bystream(aap, subjInd, spmSession, 'listspikes'));
                
                % Scans with large movement and image intensity fluctuations
                regrScans = union(TSspikes(:,1), Mspikes(:,1));
                
                % Create a delta for each of these scans
                spikes = zeros(nScans, length(regrScans));
                spikeNames = {};
                for s=1:length(regrScans),
                    spikes(regrScans(s), s) = 1;
                    spikeNames{s} = sprintf('SpikeMov%d', s);
                end;
                
                aas_log(aap, false, sprintf('%d spike regressors added for this session.', size(spikes, 2)));
                
                C = [C spikes];
                Cnames = [Cnames spikeNames];
                
            end
            
            % Remove scans if needed
            if ~isempty(settings.numDummies) && settings.numDummies > 0
                fprintf('Removing 1st %d scans.\n', settings.numDummies);
                C = C(settings.numDummies+1:end, :);
                thisSessFiles = thisSessFiles(settings.numDummies+1:end, :);
                nScans = size(thisSessFiles, 1);
            end
            
            SPM.nscan(spmSession) = nScans;
            allFiles = strvcat(allFiles, thisSessFiles);
            
            % Scale confoundes using spm_en?
            C = spm_en(C, 0); % scale confounds by sum of squares
            
            % Dimension reduction on confounds
            if svdThresh < 1
                [U,S] = spm_svd(C,0);
                if issparse(S)
                    S = full(S);
                end
                S = diag(S).^2; S = cumsum(S)/sum(S);
                Np = find(S > svdThresh); Np = Np(1);
                aas_log(aap, 0, sprintf('%d modes from %d confounds (%.2f percent variance of correlation explained)\n', Np, size(C,2), 100*S(Np)))
                X0 = full(U(:,1:Np));
                Nc = Np;
                
                % Replace confound matrix and names
                C = X0;
                Cnames = {};
                for svdReg = 1:size(X0,2)
                    Cnames{end+1} = sprintf('confound SVD %d', svdReg);
                end
            end
     
            % Check that we don't have more covariates than scans!
            if size(C,2) > SPM.nscan(spmSession)
                msg = sprintf('WARNING: More covariates (%d) than scans (%d)!', size(C,2), SPM.nscan(spmSession));
                if svdThresh==1
                    msg = [msg ' Try SVD reduction (in xml file, set svdthresh to .99 instead of 1)?'];
                end
                aas_log(aap, false, msg);
            end

            SPM.Sess(spmSession).C.C = [SPM.Sess(spmSession).C.C C];
            SPM.Sess(spmSession).C.name = Cnames;

            spmSession = spmSession + 1; % incremenet for next time
        end % looping through sessions
        
        SPM.xY.P = allFiles;
        SPM.xX.iG = nuisanceCols;
        SPM.xX.iC = interestCols;
        SPM = spm_fmri_spm_ui(SPM);
        
        % Explicit masking, if requested
        % (NB spm_fmri_spm_ui seems to reset xM.VM, so adding it here.)
        if explicitMask
            maskImg = aas_getfiles_bystream(aap, subjInd, 'native_brainmask')
            SPM.xM.VM = spm_vol(maskImg);
        end % dealing with explicit mask

        
        SPMpath = fullfile(analysisDir, 'SPM.mat');
        save(SPMpath, 'SPM');
        
        % describe outputs
        aap = aas_desc_outputs(aap, subjInd, 'firstlevel_spm', SPMpath);
        
        cd(startingDir); % go back where we started
        
        % (estimate with aamod_firstlevel_modelestimate)
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end