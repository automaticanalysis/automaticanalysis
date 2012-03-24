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
%  1 = trim second _CBUxxxxxx off (use megname instead if it exists)
%  2 = use S01, S02, S03...
%  3 = use acq_details.subjects as is
%
% examples
%   aas_getsubjpath(aap,1)   - path to first subject
%   aas_getsubjpath(aap,1,'s3')  - s3 path to first subject
%   aas_getsubjpath(aap,1,3) - path to third stage of first subject
%   aas_getsubjpath(aap,1,'s3',3) - path to third stage of first subj on s3

function [subjpath]=aas_getsubjpath(aap,i,varargin)

root=aas_getstudypath(aap,varargin{:});

switch (aap.directory_conventions.subject_directory_format)
    case 1
        if (isfield(aap.acq_details.subjects(i),'megname') && ~isempty(aap.acq_details.subjects(i).megname))
            tmp=aap.acq_details.subjects(i).megname;
        else
            numpos=findstr('CBU',aap.acq_details.subjects(i).mriname);
            tmp=aap.acq_details.subjects(i).mriname(numpos(1):length(aap.acq_details.subjects(i).mriname));
            tmp=strtok(tmp,' /\\_,.');
        end
        subjpath=fullfile(root, tmp);
    case 2
        subjpath=fullfile(root, sprintf('S%02d',i));
    case 0
        subjpath=fullfile(root, aap.directory_conventions.subject_directory_names{i});
    case 3
        subjpath=fullfile(root, aap.acq_details.subjects(i).mriname);
    otherwise
        aas_log(1,sprintf('Unknown subject directory format (aap.directory_conventions.subject_directory_format=%d',aap.directory_conventions.subject_directory_format));
end;
