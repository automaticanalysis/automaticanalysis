% Automatic analysis - get subject path
% Returns path to the data directory for the subject's analyzed data
% The particular form of the subject directories is controlled by
% aap.directory_conventions.subject_directory_format
%
% i specifies the subject number
%
% formats are:
%  0 = use directory_conventions.subject_directory_names (trimmed version produced by
%   evaluate subject names
%  1 = from data by trimming mriname or megname
%  2 = use S01, S02, S03...
%  3 = manual
%
% examples
%   aas_getsubjpath(aap,1)   - path to first subject
%   aas_getsubjpath(aap,1,'s3')  - s3 path to first subject
%   aas_getsubjpath(aap,1,3) - path to third stage of first subject
%   aas_getsubjpath(aap,1,'s3',3) - path to third stage of first subj on s3

function [subjpath]=aas_getsubjpath(aap,i,varargin)

subjpath=fullfile(aas_getstudypath(aap,varargin{:}),aap.acq_details.subjects(i).subjname);