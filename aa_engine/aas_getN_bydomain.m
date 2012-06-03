% Automatic analysis - this file determines the names of each of the
% directory levels
%   domain='subject','session' etc
%   index= number of item

function [N]=aas_getN_bydomain(aap,domain,indices)

switch (domain)
    case 'searchlight'
        Nfn=fullfile(aas_getpath_bydomain(aap,domain,[indices 1]),'searchlightN.txt');
        if exist(Nfn,'file')
            N=load(Nfn);
        else
            N=0
        end;
        
    case 'session'
        N=length(aap.acq_details.sessions);
               
    case 'subject'
        N=length(aap.acq_details.subjects);
        
end;