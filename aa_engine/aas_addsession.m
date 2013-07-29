% Automatic analysis - add session to a planned analysis
% Usually in your user script.
% Using this function supercedes separately filling 
% aap_acq_details.sessions and aap.acq_details_brukersessionnums. It is
% more convenient when you have many subjects as the correspondence between
% subject name and series numbers is more transparent.
%
% name= subject filename (may include UNIX wildcards, e.g., CBU060500*/*)
% seriesnumbers=series numbers of EPIs for this subject
% ignoreseries parameter=series numbers of any series to be ignored in the
% analysis (e.g. a repeated structural) [added by djm 20/3/06]
% specialseries= special series to be converted
%
% name= session name (e.g. Visual)

function [aap]=aas_addsession(aap,name)

% Blank template for a session entry
thissess.name=name;

% And put into acq_details, replacing a single blank entry if it exists
if (numel(aap.acq_details.sessions)==1 && isempty(aap.acq_details.sessions.name))
    aap.acq_details.sessions=thissess;
else
    doAdd = true;
    for iSess = 1:numel(aap.acq_details.sessions)
        if strcmp(aap.acq_details.sessions(iSess).name,thissess.name)
            doAdd = false;
        end
    end
    if doAdd, aap.acq_details.sessions(end+1) = thissess; end
end;