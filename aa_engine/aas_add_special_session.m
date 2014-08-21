% Add special session. Suffixes (i.e. substring after "_") will be treated
% as subsessions (similarly to the two subsessions of the fieldmap) and
% copied into subfolders. E.g.:
%
%   aap.acq_details.special_sessions.name = {'MT_baseline' 'MT_MT'}
%   aap.acq_details.subjects.specialseries = [14 15]
%
%   Series 14 and 15 will go to subjpath/MT/baseline and subjpath/MT/MT,
%   respectively.

function [aap]=aas_add_special_session(aap,name)

% Blank template for a session entry
thissess.name=name;

% And put into acq_details, replacing a single blank entry if it exists
if (length(aap.acq_details.special_sessions)==1 && isempty(aap.acq_details.special_sessions.name))
    aap.acq_details.special_sessions=thissess;
else
    aap.acq_details.special_sessions(end+1)=thissess;
end;