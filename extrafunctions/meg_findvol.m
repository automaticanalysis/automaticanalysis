% volunteer locator
%   fp - return fullpath
%
% 11 --> meg11_0011_cc710154/110210
% Tibor Auer MRC CBU Cambridge 2012-2013

function strSubj = meg_findvol(aap,subjpath,varargin)


fp = logical(cell_index(varargin,'fp'));
fdate = '';
for a = varargin
    if ~isnan(str2double(a{1})), fdate = a{1}; end
end


%% Changed to comma separated list format [RC]

% Parse comma separated list
SEARCHPATH = textscan(aap.directory_conventions.rawmegdatadir,'%s','delimiter', ':');
SEARCHPATH = SEARCHPATH{1};

% get subjname
if ~isempty(regexp(aap.directory_conventions.megsubjectoutputformat,'%s', 'once')) % string input expected
    if ~ischar(subjpath)
        aas_log(aap,true,'Second input must be a string. Check aap.directory_conventions.megsubjectoutputformat');
    end
else  % numeric input expected
    if ~isnumeric(subjpath)
        aas_log(aap,true,'Second input must be an integer. Check aap.directory_conventions.megsubjectoutputformat');
    end
end
subjpath = sprintf(aap.directory_conventions.megsubjectoutputformat,subjpath);

isFound = false;
for i = 1:numel(SEARCHPATH)
    if ~isempty(dir(fullfile(SEARCHPATH{i},subjpath)))
        isFound = true;
        break;
    end
end

if ~isFound
    aas_log(aap,false,sprintf('Subject %s* not found',subjpath));
    strSubj = '';
    return;
end

if exist(fullfile(SEARCHPATH{i},subjpath),'dir') % exact match
    strSubj = subjpath;
else % pattern
    strSubj = dir(fullfile(SEARCHPATH{i},subjpath)); strSubj = strSubj(end).name; % in case of multiple entries
end
strSubjDir = dir(fullfile(SEARCHPATH{i},strSubj));
ds = {strSubjDir(3:end).name}; 
if isempty(fdate), dmatch = ds{1}; % first entry if no date specified
else % select match
    dmatch = '';
    dind = cell_index(ds,fdate);
    if any(dind), dmatch = ds{dind(1)}; % first match
    else
        aas_log(aap,false,sprintf('Subject %s* not found',subjpath));
        strSubj = '';
        return;
    end
end

if ~isdir(fullfile(SEARCHPATH{i},strSubj,dmatch)) || isnan(str2double(dmatch)), dmatch = ''; end % test if it is a 'date_time' folder

if fp
    strSubj = fullfile(SEARCHPATH{i},strSubj,dmatch);
else
    strSubj = fullfile(strSubj,dmatch);
end
