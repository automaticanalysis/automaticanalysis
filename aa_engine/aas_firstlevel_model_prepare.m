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
    
    % Can we build an SPM model without any condition regressors?  Usually no,
    % but yes if we are using SPM for filtering, etc.
    if isfield(aap.tasklist.currenttask.settings, 'allowemptymodel') && aap.tasklist.currenttask.settings.allowemptymodel
        allowEmptyModel = 1;
    else
        allowEmptyModel = 0;
    end
    
    if ~allowEmptyModel && ( (isempty(modelnum)) && (isempty(modelCnum)) )
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

%% retrieve TR from DICOM header
% TR is manually specified (not recommended as source of error)
if isfield(aap.tasklist.currenttask.settings,'TR') && ...
        ~isempty(aap.tasklist.currenttask.settings.TR)
    SPM.xY.RT =aap.tasklist.currenttask.settings.TR;
else
    % Get TR from DICOM header checking they're the same for all sessions
    for sess=aap.acq_details.selected_sessions
        try
            DICOMHEADERS=load(aas_getfiles_bystream(aap,subj,sess,'epi_dicom_header'));
        catch
            DICOMHEADERS=load(aas_getfiles_bystream(aap,subj,sess,'epi_header')); % For backwards compatibility
        end
        try
            TR=DICOMHEADERS.DICOMHEADERS{1}.volumeTR;
        catch
        end
        % [AVG] This is for backwards compatibility!
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

%% Set up basis functions

% Start with defaults.  Some, none, or all of these might be changed in the
% .xml or the user script.
SPM.xBF.T          = 17;                % number of time bins per scan
SPM.xBF.UNITS      = 'scans';           % OPTIONS: 'scans'|'secs' for onsets
SPM.xBF.Volterra   = 1;                 % OPTIONS: 1|2 = order of convolution
SPM.xBF.name       = 'hrf';
SPM.xBF.length     = 32;                % length in seconds
SPM.xBF.order      = 1;                 % order of basis set
SPM.xBF.bf         = [];                % Custom basis functions?
% SPM.xBF.T0 is dealt with a bit further down

% Collect values from the .xml or user script
if isfield(aap.tasklist.currenttask.settings,'xBF')
    fields = fieldnames(aap.tasklist.currenttask.settings.xBF);
    for f = 1:numel(fields)
        % Overwrite the default if something is specified
        if ~isempty(aap.tasklist.currenttask.settings.xBF.(fields{f}))
            SPM.xBF.(fields{f}) = aap.tasklist.currenttask.settings.xBF.(fields{f});
        end
    end
end

SPM.xBF.dt = SPM.xY.RT / SPM.xBF.T; % Time bin length in secs

% If no custom bf is specified, remove this field so SPM uses default behaviour
if isempty(SPM.xBF.bf), SPM.xBF = rmfield(SPM.xBF, 'bf'); end

%% Allow specifying UNITS
% This should probably disappear, because UNITS is specified in xBF
if isfield(aap.tasklist.currenttask.settings,'UNITS') && ...
        ~isempty(aap.tasklist.currenttask.settings.UNITS)
    SPM.xBF.UNITS =aap.tasklist.currenttask.settings.UNITS;
    warning('UNITS should be specified in xBF.UNITS, not as it''s own setting. This setting may be removed in the future.');
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

try SPM.xGX.iGXcalc = aap.tasklist.currenttask.settings.globalscaling; catch
    SPM.xGX.iGXcalc = 'None';
end
try SPM.xVi.form = aap.tasklist.currenttask.settings.autocorrelation; catch
    SPM.xVi.form = 'AR(1)';
end

%% Adjust time bin T0 according to reference slice & slice order
%  implements email to CBU from Rik Henson 27/06/07
%  assumes timings are relative to beginning of scans
if ~isfield(SPM.xBF,'T0') || isempty(SPM.xBF.T0) % Allow T0 override in .xml, or settings
    if (usesliceorder)
        refwhen=(find(sliceorder==refslice))/(length(sliceorder));
        SPM.xBF.T = numel(sliceorder);
    else
        % Otherwise, default to halfway through the volume
        aas_log(aap,false,'No stream sliceorder found, defaulting timing to SPM.xBF.T0 to halfway through a volume.');
        refwhen = 0.5;
    end
    SPM.xBF.T0 = round(SPM.xBF.T*refwhen);
end

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

SPM.swd = anadir;