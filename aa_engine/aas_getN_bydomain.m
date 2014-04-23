% Automatic analysis - this file determines the number of parts to a given
% domain - e.g., the number of subjects at the 'subject' level, or the
% number of N-fold validations for cross-validation. 
%   domain='subject','session' etc
%   index= number of item

function [N]=aas_getN_bydomain(aap,domain,indices)

switch (domain)
    case {'searchlight_package','hyperalignment_searchlight_package'}
        N=aap.options.searchlight.Npackage;
       
    case 'session'
        N=length(aap.acq_details.sessions);

    case {'splitsession_cv_fold','splitsession_cv_fold_hyper'}
        N=aap.options.splitsession_cv.N;
               
    case {'subject','hyperalignment_subject'}
        N=length(aap.acq_details.subjects);
    case 'study'
        N=1;
        
    case {'diffusion_session', 'diffusion_session_bedpostx'}
        N=length(aap.acq_details.diffusion_sessions);
        if isempty(aap.acq_details.diffusion_sessions(1).name)
            N=N-1;
        end
        
    case 'diffusion_session_probtrackx'
        N=aap.options.probtrackx.nsplits;  
        
    case 'scan'
        N=aap.options.realtime.nscans;
end;    