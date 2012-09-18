% Automatic analysis - this file determines the names of each of the
% directory levels
%   domain='subject','session' etc
%   index= number of item

function [N]=aas_getN_bydomain(aap,domain,indices)

switch (domain)
    case {'searchlight_package','hyperalignment_searchlight_package'}
        N=aap.options.searchlight.Npackage;
       
    case 'session'
        N=length(aap.acq_details.sessions);
               
    case {'subject','hyperalignment_subject'}
        N=length(aap.acq_details.subjects);
    case 'study'
        N=0;
end;