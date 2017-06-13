% volunteer locator
%   fp - return fullpath
%
% 90952 --> CBU090952_MR09032/20090828_131456
% Tibor Auer MRC CBU Cambridge 2012-2013
%
% CHANGE HISTORY
%
% 05/17 [MSJ] added a try/catch to handle empty directory and OS X .DS_Store weirdness

function strSubj = mri_findvol(aap,subjpath,fp)

if nargin < 3, fp = false; end

%% Changed to comma separated list format [RC]

if isstruct(aap.directory_conventions.rawdatadir)
    aas_log(aap,true,'Structure format for rawdatadir no longer supported - use comma separated list');
end;

% Parse comma separated list

SEARCHPATH = textscan(aap.directory_conventions.rawdatadir,'%s','delimiter', ':');
SEARCHPATH = SEARCHPATH{1};

% get subjname
if ~isempty(regexp(aap.directory_conventions.subjectoutputformat,'%s', 'once')) % string input expected
	if ~ischar(subjpath)
		aas_log(aap,true,'Second input must be a string. Check aap.directory_conventions.subjectoutputformat');
	end
else  % numeric input expected
	if ~isnumeric(subjpath)
    	aas_log(aap,true,'Second input must be an integer. Check aap.directory_conventions.subjectoutputformat');
	end
end
subjpath = sprintf(aap.directory_conventions.subjectoutputformat,subjpath);

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

strSubjDir = dir(fullfile(SEARCHPATH{i},strSubj));

% test if f1 is a 'date_time' folder 

% numel(strSubjDir) will only be 2 ("." and "..") if the directory is empty. 
% Also, strSubjDir(3) might be ".DS_Store" if OS X (it might not be -- it 
% depends on how the folder was created). 

try
	f1 = strSubjDir(3).name;	% first directory entry?
catch
	f1 = '';
end
if (strcmp(f1,'.DS_Store'))
	try
		f1 = strSubjDir(4).name;
	catch
		f1 = '';
	end
end
try
    junk = datenum(f1,'yyyymmdd_HHMMSS'); 
catch
    f1 = '';
end

if fp
    strSubj = fullfile(SEARCHPATH{i},strSubj,f1);
else
    strSubj = fullfile(strSubj,f1);
end
