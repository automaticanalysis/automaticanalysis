% Automatic analysis - get session path
% Returns path to the data directory for a particular session of a particular subject's analyzed data
% The particular form of the subject directories is controlled by
% aap.directory_conventions.subject_directory_format
%
% examples:
%  pth=aas_getsesspath(aap,1,2);  % subject 1, session 2
%  pth=aas_getstudypath(aap,1,2,'s3');   % path on s3 
%  pth=aas_getstudypath(aap,1,2,'s3',4);   % path on s3 for module 4

function [sesspath]=aas_getsesspath(aap,i,j,varargin)

global defaults
if ~isstruct(defaults) || ~isfield(defaults,'modality')
    aas_log(aap,0,'WARNING:defaults.modality is not set; FMRI is assumed');
    defaults.modality = 'FMRI'; % default modality
end

switch defaults.modality
    case 'FMRI'
        sesspath=fullfile(aas_getsubjpath(aap,i,varargin{:}),aap.acq_details.sessions(j).name);
    case 'EEG'
        sesspath=fullfile(aas_getsubjpath(aap,i,varargin{:}),aap.acq_details.meg_sessions(j).name);
end