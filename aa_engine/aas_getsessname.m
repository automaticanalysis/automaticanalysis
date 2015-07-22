function [nme]=aas_getsessname(aap,i,j)

if ~isfield(aap.spm.defaults,'modality')
    aas_log(aap,0,'WARNING:modality is not set; (F)MRI is assumed');
    aap.spm.defaults.modality = 'FMRI'; % default modality
end

switch aap.spm.defaults.modality
    case 'FMRI'
        sessions = aap.acq_details.sessions;
    case 'EEG'
        sessions = aap.acq_details.meg_sessions;
end

if ~isempty(sessions(j).name)
    nme = sessions(j).name;
else
    nme='(unknown)';
end;

nme=[aas_getsubjname(aap,i),'; session ',nme];
