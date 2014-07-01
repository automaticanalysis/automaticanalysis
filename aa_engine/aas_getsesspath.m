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

sesspath=fullfile(aas_getsubjpath(aap,i,varargin{:}),aap.acq_details.sessions(j).name);