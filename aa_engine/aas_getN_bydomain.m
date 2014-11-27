% Automatic analysis - this file determines the number of parts to a given
% domain - e.g., the number of subjects at the 'subject' level, or the
% number of N-fold validations for cross-validation.
%
% Updated to allow different N, depending on the indices that are passed
% in. Useful for cases where N is variable; e.g., if a subject is missing a
% session.
%
% Now we also return the indices of those parts. For example, if there are
% three sessions per subject (N=3, I=[1 2 3]), but one subject is missing
% the middle session (N=2, for that subject), then the session indices for
% that subject should be [1 3]
%
%   domain = 'subject','session' etc
%   index = number of item
% -----------------------------------------------------------------------
% Updated by cwild 2014-05-12 to allow variable N for sessions. 

function [N, I]=aas_getN_bydomain(aap,domain,indices)

switch (domain)
    case {'searchlight_package','hyperalignment_searchlight_package'}
        N=aap.options.searchlight.Npackage;
        I=1:N;
        
    case {'session','meg_session'}
        
        switch domain
            case 'session'
                sessions = aap.acq_details.sessions;
                seriesnumbers = aap.acq_details.subjects(indices(1)).seriesnumbers;
            case 'meg_session'
                sessions = aap.acq_details.meg_sessions;
                seriesnumbers = aap.acq_details.subjects(indices(1)).megseriesnumbers;
        end
                
        if isempty(indices)
            N = length(sessions);
            I = 1 : N;
        else
            if iscell(seriesnumbers)
                N = cellfun(@(x) ~isempty(x), seriesnumbers);
                I = find(N);
                N = sum(N);
            else
                N=sum(seriesnumbers > 0);
                I=find(seriesnumbers > 0);
            end
        end
        
    case {'splitsession_cv_fold','splitsession_cv_fold_hyper'}
        N=aap.options.splitsession_cv.N;
        I=1:N;
               
    case {'subject','hyperalignment_subject'}
        N=length(aap.acq_details.subjects);
        I=1:N;
        
    case 'study'
        N=1;
        I=1:N;
        
    case {'diffusion_session', 'diffusion_session_bedpostx'}
        N=length(aap.acq_details.diffusion_sessions);
        if isempty(aap.acq_details.diffusion_sessions(1).name)
            N=N-1;
        end
        I=1:N;
        
    case 'diffusion_session_probtrackx'
        N=aap.options.probtrackx.nsplits;  
        I=1:N;
        
    case 'scan'
        N=aap.options.realtime.nscans;
        I=1:N;
end;    