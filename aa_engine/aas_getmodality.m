function modality = aas_getmodality(aap)
modality = '';
if isfield(aap.tasklist.currenttask,'modality'), modality = aap.tasklist.currenttask.modality; end
if isempty(modality) && isfield(aap.spm.defaults,'modality'), modality = aap.spm.defaults.modality; end

% MRI (old) is meaningless
if strcmp(modality,'MRI'), modality = ''; end

% last resort --> try modulename
if isempty(modality)
    if strfind(aap.tasklist.currenttask.name,'_epi'), modality = 'FMRI'; end
    if strfind(aap.tasklist.currenttask.name,'_diffusion'), modality = 'DWI'; end
    if strfind(aap.tasklist.currenttask.name,'_MTI'), modality = 'MTI'; end
    if strfind(aap.tasklist.currenttask.name,'_meg'), modality = 'MEG'; end
end

% default
if isempty(modality)
%     aas_log(aap,0,'WARNING:modality cannot be determined; (F)MRI is assumed');
    modality = 'FMRI'; % default modality
end

end