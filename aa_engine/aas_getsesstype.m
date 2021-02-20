function session = aas_getsesstype(aap)

%% Modality
modality = '';
if isfield(aap.tasklist.currenttask,'modality'), modality = aap.tasklist.currenttask.modality; end
if isempty(modality) && isfield(aap.spm.defaults,'modality'), modality = aap.spm.defaults.modality; end

% MRI (old) is meaningless
if strcmp(modality,'MRI'), modality = ''; end

% if module has specific session domain
switch aap.tasklist.currenttask.domain
    case 'diffusion_session'
        modality = 'DWI';
    case 'meeg_session'
        modality = 'MEEG';
    case 'special_session'
        modality = 'X';
    otherwise
        % ignore generic session
end

% last resort --> try modulename
if isempty(modality)
    if strfind(aap.tasklist.currenttask.name,'_epi'), modality = 'FMRI'; end
    if strfind(aap.tasklist.currenttask.name,'_diffusion'), modality = 'DWI'; end
    if strfind(aap.tasklist.currenttask.name,'_MTI'), modality = 'X'; end
    if strfind(aap.tasklist.currenttask.name,'_ASL'), modality = 'X'; end
    if strfind(aap.tasklist.currenttask.name,'_meeg'), modality = 'MEEG'; end
    if strfind(aap.tasklist.currenttask.name,'_meg'), modality = 'MEG'; end
    if strfind(aap.tasklist.currenttask.name,'_eeg'), modality = 'EEG'; end
end

% default
if isempty(modality)
%     aas_log(aap,0,'WARNING:modality cannot be determined; (F)MRI is assumed');
    modality = 'FMRI'; % default modality
end

%% Session
switch modality
    case 'FMRI'
        session = 'session';
    case 'DWI'
        session = 'diffusion_session';
    case {'MEEG' 'MEG' 'EEG'}
        session = 'meeg_session';
    case 'X'
        session = 'special_session';
end
end