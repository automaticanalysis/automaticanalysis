function modality = aas_getmodality(aap)
modality = '';
if isfield(aap.tasklist.currenttask,'modality'), modality = aap.tasklist.currenttask.modality; end
if isempty(modality) && isfield(aap.spm.defaults,'modality'), modality = aap.spm.defaults.modality; end

if isempty(modality)
    aas_log(aap,0,'WARNING:modality is not set; (F)MRI is assumed');
    modality = 'FMRI'; % default modality
end

% backward compatibility
if strcmp(modality,'MRI'), modality = 'FMRI'; end
end