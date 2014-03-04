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
% Modified for aa by Rhodri Cusack MRC CBU Mar 2006Aug 2007
% Thanks to Rik Henson for various suggestions

function [aap,resp]=aamod_firstlevel_model_MVPaa(aap,task,subj)

resp='';

switch task
    case 'report'
        
    case 'doit'
        % Get subject directory
        cwd=pwd;
        
        % Prepare basic SPM model...
        [SPM, anadir, files, allfiles, model, modelC] = aas_firstlevel_model_prepare(aap, subj);
        
        % Get all the nuisance regressors...
        [movementRegs, compartmentRegs, physiologicalRegs, spikeRegs, GLMDNregs] = ...
            aas_firstlevel_model_nuisance(aap, subj, files);
        
        %% Set up CORE model
        coreSPM = SPM;
        cols_nuisance=[];
        cols_interest=[];
        currcol=1;
        
        sessnuminspm=0;
        
        for sess = aap.acq_details.selected_sessions
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
        
        SPM.xY.P = allfiles;
        SPMdes = spm_fmri_spm_ui(SPM);
        
        %% DIAGNOSTIC
        mriname = aas_prepare_diagnostic(aap, subj);
        try
            saveas(1, fullfile(aap.acq_details.root, 'diagnostics', ...
                [mfilename '__' mriname '.fig']));
        catch
        end
        
        % now check real covariates and nuisance variables are
        % specified correctly
        SPMdes.xX.iG=cols_nuisance;
        SPMdes.xX.iC=cols_interest;
        
        % Turn off masking if requested
        if ~aap.tasklist.currenttask.settings.firstlevelmasking
            SPMdes.xM.I=0;
            SPMdes.xM.TH=-inf(size(SPMdes.xM.TH));
        end
        
        spm_unlink(fullfile('.', 'mask.img')); % avoid overwrite dialog
        set(0,'RecursionLimit',2000);
        SPMest = spm_spm(SPMdes);
        
        %% REDO MODEL WITH Mumford/Poldrak method...
        
        % Find out how large the model{sess} should be (per session)
        eventNumber = [];
        sessNumber = [];
        
        for sess = aap.acq_details.selected_sessions
            eventNumber = [eventNumber 1:size(model{sess}.event,2)];
            sessNumber = [sessNumber sess*ones(1,size(model{sess}.event,2))];
        end
        sessRegs = 1:max(eventNumber);
        
        % Loop over regressors and do Mumford/Poldrack modelling
        
        
        switch aap.tasklist.currenttask.settings.parallel
            case 'torque'
                for numReg = sessRegs
                    aapCell{numReg}                = aap;
                    anadirCell{numReg}             = anadir;
                    coreSPMCell{numReg}            = coreSPM;
                    filesCell{numReg}              = files;
                    allfilesCell{numReg}           = allfiles;
                    
                    modelCell{numReg}              = model;
                    modelCCell{numReg}             = modelC;
                    eventNumberCell{numReg}        = eventNumber;
                    sessNumberCell{numReg}         = sessNumber;
                    numRegCell{numReg}             = numReg;
                    nDestCell{numReg}              = SPMest.xX.iC(eventNumber==numReg);
                    
                    movementRegsCell{numReg}       = movementRegs;
                    compartmentRegsCell{numReg}    = compartmentRegs;
                    physiologicalRegsCell{numReg}  = physiologicalRegs;
                    spikeRegsCell{numReg}     = spikeRegs;
                end
                
                qsubcellfun(@aas_firstlevel_model_mumford, ...
                    aapCell, anadirCell, coreSPMCell, filesCell, allfilesCell, ...
                    modelCell, modelCCell, eventNumberCell, sessNumberCell, numRegCell, nDestCell, ...
                    movementRegsCell, compartmentRegsCell, physiologicalRegsCell, spikeRegsCell, ...
                    'memreq', 3.7 * (1024^3), ... % NOT DYNAMIC YET!!!
                    'timreq', 1 * (3600), ...
                    'stack', 1 ...
                    );
            otherwise % serial/none
                for numReg = sessRegs
                    nDest = SPMest.xX.iC(eventNumber==numReg);
                    aas_firstlevel_model_mumford(aap, anadir, coreSPM, files, allfiles, ...
                        model, modelC, eventNumber, sessNumber, numReg, nDest, ...
                        movementRegs, compartmentRegs, physiologicalRegs, spikeRegs)
                end
        end
        
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
        if ~aap.tasklist.currenttask.settings.firstlevelmasking
            otherfiles={'ResMS.hdr','ResMS.img','RPV.hdr','RPV.img'};
        else
            otherfiles={'mask.hdr','mask.img','ResMS.hdr','ResMS.img','RPV.hdr','RPV.img'};
        end
        for otherind=1:length(otherfiles)
            betafns=strvcat(betafns,fullfile(anadir,otherfiles{otherind}));
        end
        aap=aas_desc_outputs(aap,subj,'firstlevel_betas',betafns);
        
        %% DIAGNOSTICS...
        firstlevelmodelStats(anadir, [], fullfile(anadir, 'mask.img'))
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
