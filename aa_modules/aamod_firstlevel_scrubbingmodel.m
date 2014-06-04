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
         
        % Prepare basic SPM model...
        [SPM, anadir, files, allfiles, model, modelC] = aas_firstlevel_model_prepare(aap, subjInd);
        cd(SPM.swd);
        
        TR = SPM.xY.RT;
        
        % Loop through sessions and specify model
        nuisanceCols = [];
        interestCols = [];
        currentCol = 1;
        
        spmSession = 1;
        
        allFiles = [];
 
        for thisSess = aap.acq_details.selected_sessions       
            
            % Set up model, primarily for experimental effects. We are
            % going to add the nuisance columns on our own...
            [SPM, interestCols, nuisanceCols, currentCol] = ...
                aas_firstlevel_model_define(aap, thisSess, spmSession, SPM, model, modelC, ...
                                                             interestCols, nuisanceCols, currentCol, ...
                                                             [], [], [], [], []);

            nScans = size(files{thisSess},1);
   
            % Set HP filter, if there is one
            SPM.xX.K(spmSession).HParam = highpassFilter; 

            % Confounds/covariates of no interest/whatever you want to call
            % them. These are collected together in case one wants to do
            % dimensionality reduction.
            
            C = []; % confounds
            Cnames = []; % names

            % Only get compFiles if requested (allows more flexibility in
            % XML files)
            
            if settings.includeCSF || settings.includeWM
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
                if settings.includeWM
                    fprintf('Adding WM signal as covariate.\n');
                    C = [C compSignals.compTC(:,2)];
                    Cnames = [Cnames 'WM'];
                end
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
            
%             settings.includeSpikes=0; % hack because includeSpikes was a cell?!? -JP
            if settings.includeSpikes && any(strcmp('listspikes', aap.tasklist.currenttask.inputstreams.stream))
                
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
                files{thisSess} = files{thisSess}(settings.numDummies+1:end, :);
                nScans = size(files{thisSess}, 1); 
            end
            SPM.nscan(thisSess) = nScans;
            
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
            
            nuisanceCols = [nuisanceCols [1:size(C,2)]+(currentCol-1)];
            currentCol = currentCol + size(C, 2);
            
            SPM.Sess(spmSession).C.C = [SPM.Sess(spmSession).C.C C];
            SPM.Sess(spmSession).C.name = Cnames;

            spmSession = spmSession + 1;
        end % looping through sessions
        
        % Add the images
        SPM.xY.P = char(files{:});
        
        if settings.nuisanceconditions
            nuisanceCols = [nuisanceCols interestCols];
            interestCols = [];
        end
        
        SPM.xX.iG = nuisanceCols;
        SPM.xX.iC = interestCols;
        
        SPM = spm_fmri_spm_ui(SPM);
        
        % Explicit masking, if requested
        % (NB spm_fmri_spm_ui seems to reset xM.VM, so adding it here.)
        if explicitMask
            maskImg = aas_getfiles_bystream(aap, subjInd, 'native_brainmask')
            SPM.xM.VM = spm_vol(maskImg);
        end % dealing with explicit mask

        SPMpath = fullfile(anadir, 'SPM.mat');
        save(SPMpath, 'SPM');
        
        % describe outputs
        aap = aas_desc_outputs(aap, subjInd, 'firstlevel_spm', SPMpath);
        
        cd(startingDir); % go back where we started

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end