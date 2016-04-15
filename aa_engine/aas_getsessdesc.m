function [nme]=aas_getsessdesc(aap,i,j)

switch aas_getmodality(aap)
    case 'FMRI'
        sessions = aap.acq_details.sessions;
    case 'DWI'
        sessions = aap.acq_details.diffusion_sessions;
    case {'MEG' 'EEG'}
        sessions = aap.acq_details.meg_sessions;
    case {'MTI'}
        sessions = aap.acq_details.special_sessions;
end

if ~isempty(sessions(j).name)
    nme = sessions(j).name;
else
    nme='(unknown)';
end;

nme=[aas_getsubjdesc(aap,i),'; session ',nme];
