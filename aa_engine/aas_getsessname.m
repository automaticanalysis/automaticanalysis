function [nme]=aas_getsessname(aap,i,j)

switch aas_getmodality(aap)
    case 'FMRI'
        sessions = aap.acq_details.sessions;
    case 'DWI'
        sessions = aap.acq_details.diffusion_sessions;
    case {'MEG' 'EEG'}
        sessions = aap.acq_details.meg_sessions;
end

if ~isempty(sessions(j).name)
    nme = sessions(j).name;
else
    nme='(unknown)';
end;

nme=[aas_getsubjname(aap,i),'; session ',nme];
