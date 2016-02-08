% Automatic analysis - this file determines the names of each of the
% directory levels
%   domain='subject','session' etc
%   index= number of item

function [directory]=aas_getdirectory_bydomain(aap,domain,index)

switch (domain)
    case 'searchlight_package'
        directory=sprintf('searchlight_package_%d',index);
    
    case 'hyperalignment_searchlight_package'
        directory=sprintf(['hyperalignment_searchlight_packages' filesep '%d'],index);
        
    case {'splitsession_cv_fold','splitsession_cv_fold_hyper'}
        directory=sprintf('splitsession_cv_fold_%d',index);

    case 'diffusion_session_probtrackx'
        directory=sprintf('probtrackx_%d',index);

    case 'session'
        directory=aap.acq_details.sessions(index).name;
        
    case 'isc_session'
        directory = ['isc_' aap.acq_details.sessions(index).name];

    case 'scan'
        directory=sprintf('scan_%d',index);
    
    case 'meg_session'
        directory=aap.acq_details.meg_sessions(index).name;
        
    case 'diffusion_session'
        directory=aap.acq_details.diffusion_sessions(index).name;

    case 'diffusion_session_bedpostx'
        directory=[aap.acq_details.diffusion_sessions(index).name '.bedpostX'];
        
    case {'subject','hyperalignment_subject'}
        directory=aap.acq_details.subjects(index).subjname;  
end