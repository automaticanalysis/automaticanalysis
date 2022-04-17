% Download one of the predefined datasets into aap.directory_conventions.rawdatadir.
%
% Only does the download if aap.directory_conventions.rawdatadir does not yet exist or is empty.
%
% Inputs:
%  - aap: The aap structure
%  - dataset_id: One of the following values
%    + 'aa_demo'
%    + 'ds114_test2'
%    + 'MoAEpilot'
%    + 'ds000114', N.B.: you may want to specify subset because the whole dataset is large
%    + 'ds002737', N.B.: you may want to specify subset because the whole dataset is large
%    + 'LEMON_EEG', N.B.: you may want to specify subset because the whole dataset is large
%    + 'LEMON_MRI', N.B.: you may want to specify subset because the whole dataset is large
%
% aa_downloaddemo(aap, dataset_id)
function aa_downloaddemo(aap, dataset_id, subset)

%% Inputs checking
demodir = aap.directory_conventions.rawdatadir;
% When used in aas_log messages, escape backward slashes from windows paths.
logsafe_path = strrep(demodir, '\', '\\');

% Check aap.directory_conventions.rawdatadir
sources = strsplit(demodir, pathsep);
if length(sources)>1
    % only want one rawdatadir for downloaddemo
    msg = sprintf('ERROR: For use with aa_downloaddemo, aap.directory_conventions.rawdatadir (%s) must specify exactly one directory.', logsafe_path);
    aas_log(aap, true, msg);
end

% Check dataset_id
datasets = datasetClass.empty;
for d = jsondecode(fileread('datasets.json'))'
    par = struct2cell(d);
    datasets(end+1) = datasetClass(par{:});
end
IDs = {datasets.ID};
ID_ind = strcmp(dataset_id, IDs);
if sum(ID_ind) ~= 1
    IDs_str = strjoin(IDs,', ');
    msg = sprintf('ERROR: Expected exactly one match for input dataset_id (%s) in list of known datasets %s', dataset_id, IDs_str);
    aas_log(aap, true, msg);
end
dataset = datasets(ID_ind);

%% Download if not already has data
if ~exist(fullfile(demodir),'dir') ... % Does not exist yet
        || length(dir(demodir))<3 % Directory is empty (only . and .. entries in dir listing)

    [mkdir_status, mkdir_msg] = mkdir(demodir); % Create if needed
    if ~mkdir_status
        msg = sprintf('ERROR: failed to create directory %s, due to: %s', logsafe_dir, mkdir_msg);
        aas_log(aap, true, msg);
    end

    aas_log(aap, false, ['INFO: downloading demo data to ' logsafe_path]);

    % Download and unpack the data to a temp dir first
    if nargin == 3, dataset.subset = subset; end
    dataset.download(demodir);    
else
    msg = sprintf('INFO: aa_downloaddemo: Directory %s is already non-empty, skipping data download', logsafe_path);
    aas_log(aap, false, msg);
end

end
