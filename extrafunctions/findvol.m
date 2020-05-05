% volunteer locator
%   fp - return fullpath
%
% 90952 --> CBU090952_MR09032/20090828_131456
% Tibor Auer MRC CBU Cambridge 2012-2013
%
% CHANGE HISTORY
%
% 05/17 [MSJ] added a try/catch to handle empty directory and OS X .DS_Store weirdness

function strSubj = findvol(aap,modality,subjpath,fdate,fp)

DTFORMAT = 'yyyymmdd_HHMMSS';

switch modality
    case 'mri'
        rawdatadir = aap.directory_conventions.rawdatadir;
        subjectoutputformat = aap.directory_conventions.subjectoutputformat;
    case 'meeg'
        rawdatadir = aap.directory_conventions.rawmeegdatadir;
        subjectoutputformat = aap.directory_conventions.meegsubjectoutputformat;
end

%% Changed to comma separated list format [RC]
if isstruct(rawdatadir)
    aas_log(aap,true,'Structure format for rawdatadir no longer supported - use comma separated list');
end

% Parse comma separated list
SEARCHPATH = textscan(rawdatadir,'%s','delimiter', ':');
SEARCHPATH = SEARCHPATH{1};

% get subjname
if ~isempty(regexp(subjectoutputformat,'%s', 'once')) % string input expected
	if ~ischar(subjpath)
		aas_log(aap,true,'Second input must be a string. Check aap.directory_conventions.subjectoutputformat');
	end
else  % numeric input expected
	if ~isnumeric(subjpath)
    	aas_log(aap,true,'Second input must be an integer. Check aap.directory_conventions.subjectoutputformat');
	end
end
subjpath = sprintf(subjectoutputformat,subjpath);

isFound = false;
for i = 1:numel(SEARCHPATH)
    if ~isempty(dir(fullfile(SEARCHPATH{i},subjpath)))
        isFound = true;
        break;
    end
end

if ~isFound
    aas_log(aap,true,sprintf('Subject %s* not found',subjpath));
    strSubj = '';
    return;
end

if exist(fullfile(SEARCHPATH{i},subjpath),'dir') % exact match
    strSubj = subjpath;
else % pattern
    strSubj = dir(fullfile(SEARCHPATH{i},subjpath)); strSubj = strSubj(end).name; % in case of multiple entries
end

% regexp to find files/folders that start with alphanumeric characters (ignore . files)
strSubjDir = spm_select('List',fullfile(SEARCHPATH{i},strSubj),'dir','^[a-zA-Z0-9]*');
% if there is no subdirectory
if isempty(strSubjDir), strSubjDir = spm_select('List',fullfile(SEARCHPATH{i},strSubj),'^[a-zA-Z0-9]*'); end
if isempty(strSubjDir)
    aas_log(aap,true,sprintf('nothing found for path %s',...
        fullfile(SEARCHPATH{i},strSubj)));
end

% handle CBU-style sub-directories with date formats
strSubjDir = cellstr(strSubjDir);
des = cellfun(@(x) datenum_safe(x,DTFORMAT),strSubjDir);

if ~any(des)
    if ~isempty(fdate) % subfolder of specific date required
        aas_log(aap,false,sprintf('Subject %s* on date %s not found',subjpath,fdate));
        strSubj = '';
        return;
    else % no sub-folder required
        subfolder = '';
    end
else
    if ~isempty(fdate) % subfolder of specific date required
        ind = find(des == datenum(fdate,DTFORMAT),1,'first'); % first matching subfolder (we do not even expect more)
        if isempty(ind)
            aas_log(aap,false,sprintf('Subject %s* on date %s not found',subjpath,fdate));
            strSubj = '';
            return;
        else
            subfolder = strSubjDir{ind};
        end
    else % no sub-folder required
        strSubjDir = strSubjDir(logical(des));
        subfolder = strSubjDir{1}; % first subfolder (if there is more, they are ignored)
    end
end


if fp
    strSubj = fullfile(SEARCHPATH{i},strSubj,subfolder);
else
    strSubj = fullfile(strSubj,subfolder);
end
end

function de = datenum_safe(datestr,dtformat)
try
    de = datenum(datestr,dtformat);
catch
    de = 0;
end
end
