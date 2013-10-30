function [SPM, anadir, files, allfiles, model, modelC] = ...
    aas_firstlevel_model_prepare(aap, subj)

subjname = aap.acq_details.subjects(subj).mriname;

%% Get data (EPI) files, and the models...
files = cell(size(aap.acq_details.sessions));
allfiles='';
model = cell(size(aap.acq_details.sessions));
modelC = cell(size(aap.acq_details.sessions));

for sess = aap.acq_details.selected_sessions
    files{sess} = aas_getfiles_bystream(aap,subj,sess,'epi');
    if isfield(aap.options, 'NIFTI4D') && aap.options.NIFTI4D % 4D
        V = spm_vol(files{sess});
        f0 = files{sess};
        files{sess} = '';
        for f = 1:numel(V)
            files{sess} = strvcat(files{sess},[f0 ',' num2str(V(f).n(1))]);
        end
    end
    allfiles = strvcat(allfiles,files{sess});
end

%%%%%%%%%%%%%%%%%%%
%% Modeling      %%
%%%%%%%%%%%%%%%%%%%
for sess = aap.acq_details.selected_sessions
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
    modelnum = find(sessmatches & subjmatches);
    
    %% Get modelC (covariate) data from aap
    if isfield(aap.tasklist.currenttask.settings, 'modelC')
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
        modelCnum = find(sessmatches & subjmatches);
    else
        modelCnum = [];
    end
    
    %% Check that we have at least one model of interest (normal or covariate)
    if (length(modelnum)>1) || (length(modelCnum)>1)
        aas_log(aap,true,sprintf('Error while getting model details as more than one specification for subject %s session %s',subjname,aap.acq_details.sessions(sess).name));
    end
    if (isempty(modelnum)) && (isempty(modelCnum))
        aas_log(aap,true,'Cannot find model specification. Check either user script or aamod_firstlevel_model.xml');
    end
    
    if ~isempty(modelnum)
        model{sess} = aap.tasklist.currenttask.settings.model(modelnum);
    else
        model{sess} = [];
    end
    if ~isempty(modelCnum)
        modelC{sess} = aap.tasklist.currenttask.settings.modelC(modelnum);
    else
        modelC{sess} = [];
    end
    
end

%%%%%%%%%%%%%%%%%%%
%% SPM structure %%
%%%%%%%%%%%%%%%%%%%
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

%% Allow specifying UNITS
if isfield(aap.tasklist.currenttask.settings,'UNITS') && ...
        ~isempty(aap.tasklist.currenttask.settings.UNITS)
    SPM.xBF.UNITS =aap.tasklist.currenttask.settings.UNITS;
end
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
        end
        % [AVG] This is for backwards compatibility!
        % [TA] volumeTR might be valid but empty!
        if ~exist('TR','var') || isempty(TR)
            TR=DICOMHEADERS.DICOMHEADERS{1}.RepetitionTime/1000;
        end

        if (sess==aap.acq_details.selected_sessions(1))
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
usesliceorder = aas_stream_has_contents(aap,'sliceorder');
if (usesliceorder)
    for sess=aap.acq_details.selected_sessions
        sliceorderstruct=load(aas_getfiles_bystream(aap,subj,sess,'sliceorder'));
        if (sess==aap.acq_details.selected_sessions(1))
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
    refwhen=(find(sliceorder==refslice))/(length(sliceorder));
    SPM.xBF.T = numel(sliceorder);
else
    aas_log(aap,false,'No stream sliceorder found, defaulting timing to SPM.xBF.T0=1 in model');
    refwhen=1;
end
SPM.xBF.T0 = round(SPM.xBF.T*refwhen);

%% Deal with extraparameters. Not needed any more, as
% aap.directory_conventions.stats_singlesubj
% can have module specific value, but kept for backwards compatability
if (isfield(aap.tasklist.currenttask.extraparameters,'stats_suffix'))
    stats_suffix=aap.tasklist.currenttask.extraparameters.stats_suffix;
else
    stats_suffix=[];
end

anadir = fullfile(aas_getsubjpath(aap,subj), [aap.directory_conventions.stats_singlesubj stats_suffix]);
if ~exist(anadir,'dir')
    mkdir(anadir);
end