% volunteer locator
%   fp - return fullpath
%
% 90952 --> CBU090952_MR09032/20090828_131456
% Tibor Auer MRC CBU Cambridge 2012-2013

function strSubj = mri_findvol(aap,subjpath,fp)

if nargin < 3, fp = false; end


%% Changed to comma separated list format [RC]
if isstruct(aap.directory_conventions.rawdatadir)
    aas_log(aap,true,'Structure formate for rawdatadir no longer supported - use comma separated list');
end;

% Parse comma separated list
SEARCHPATH={};
rem=aap.directory_conventions.rawdatadir;
while length(rem)>0
    [pth rem]=strtok(rem,':');
    SEARCHPATH{end+1}=pth;
end;

% get subjname
if isnumeric(subjpath)
    subjpath = sprintf(aap.directory_conventions.subjectoutputformat,subjpath);
elseif ischar(subjpath) % mriname
    subjpath = [aas_mriname2subjname(subjpath) '*'];
end

isFound = false;
for i = 1:numel(SEARCHPATH)
    if ~isempty(dir(fullfile(SEARCHPATH{i},subjpath)))
        isFound = true;
        break;
    end
end

if ~isFound
    fprintf('Subject %s* not found',subjpath);
    strSubj = '';
    return;
end

if exist(fullfile(SEARCHPATH{i},subjpath),'dir') % exact match
    strSubj = subjpath;
else % pattern
    strSubj = dir(fullfile(SEARCHPATH{i},subjpath)); strSubj = strSubj(end).name; % in case of multiple entries
end
strSubjDir = dir(fullfile(SEARCHPATH{i},strSubj));
f1 = strSubjDir(3).name; % first entry
try
    junk = datenum(f1,'yyyymmdd_HHMMSS'); % test if it is a 'date_time' folder
catch
    f1 = '';
end
if fp
    strSubj = fullfile(SEARCHPATH{i},strSubj,f1);
else
    strSubj = fullfile(strSubj,f1);
end
