% Download one of the predefined datasets into aap.directory_conventions.rawdatadir.
%
% Only does the download if aap.directory_conventions.rawdatadir does not yet exist or is empty.
%
% Inputs:
%  - aap: The aap structure
%  - dataset_id: One of the following values
%    + 'aa_demo'
%    + 'ds000114'
%
% aa_downloaddemo(aap, dataset_id)
function aa_downloaddemo(aap, dataset_id)

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
datasets = aa_downloaddemo_datasets();
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
    tgz_filename = [tempname dataset.filetype];
    tgz_filename = websave(tgz_filename, dataset.URL);
    unpack_dir = tempname;
    untar(tgz_filename, unpack_dir);

    % Depending on the contents of the downloaded tgz, different datasets
    % need different postprocessing to have the correct contents in demodir
    dataset.postprocessing(unpack_dir, demodir);

    % Check success
    if ~exist(fullfile(demodir),'dir') ... % Does not exist yet
        || length(dir(demodir))<3 % Directory is empty (only . and .. entries in dir listing)
        % assert(exist(sources{demoind},'dir')~=0);
    else
        aas_log(aap, false, 'INFO: done');
    end

    % delete the downloaded archive
    delete(tgz_filename);
else
    msg = sprintf('INFO: aa_downloaddemo: Directory %s is already non-empty, skipping data download', logsafe_path);
    aas_log(aap, false, msg);
end

end
