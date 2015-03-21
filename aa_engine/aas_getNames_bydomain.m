% Automatic analysis - this file returns the names of levels of a domain.
% For example, if domain is 'subject', we return all the subject names. If
% the domain is 'session' we return the session names. As you add new
% domains, add the appropriate thing here.
%   domain='subject','session' etc
%
% Added by CW: 2014-04-02
%
function [names] = aas_getNames_bydomain(aap, domains)

if ~iscell(domains), domains = {domains}; end

for d = 1 : length(domains)
    
    switch (domains{d})
        
        case 'session'
            names{d} = {aap.acq_details.sessions.name};

        case 'isc_session'
            for sessind=1:length(aap.acq_details.sessions)
                names{d}{sessind} = ['isc_' aap.acq_details.sessions(sessind).name];
            end;
            
        case {'subject'}
            names{d} = {aap.acq_details.subjects.mriname};
            
        case 'study'
            names{d} = {aap.directory_conventions.analysisid};
            
        case 'meg_session'
            names{d} = {aap.acq_details.meg_sessions.name};
            
        case {'diffusion_session', 'diffusion_session_bedpostx'}
            names{d} = {aap.acq_details.diffusion_sessions.name};
            
        case 'scan'
            names{d} = {};
            
        case 'searchlight_package'
            names{d} = {};
            
        case 'splitsession_cv_fold'
            names{d} = {};
            
        case 'diffusion_session_probtrackx'
            names{d} = {};      

        case 'hyperalignment_searchlight_package'
            names{d} = {};
            
        case 'hyperalignment_subject'
            names{d} = {};
            
        case 'splitsession_cv_fold_hyper'
            names{d} = {};
            
        otherwise
            aas_log(aap, 1, sprintf('Invalid domain ''%s'', perhaps NYI in this function', domains{d}));
    end
    
    if all(cellfun(@(x) isempty(x), names{d})),
        names{d} = {};
    end
    
end

% if length(names) == 1, names = names{1}; end