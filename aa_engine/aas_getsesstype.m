function session = aas_getsesstype(aap)

switch aas_getmodality(aap)
    case 'FMRI'
        session = 'session';
    case 'DWI'
        session = 'diffusion_session';
    case {'MEG' 'EEG'}
        session = 'meg_session';
    case {'MTI'}
        session = 'special_session';
end
end