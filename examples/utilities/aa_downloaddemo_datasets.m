function datasets = aa_downloaddemo_datasets()
% aa_downloaddemo_datasets returns a list of datasets that can be selected in aa_downloaddemo

% NOTE: Separate function, to enable easy testing that all listed datasets are still valid.
% NOTE: postprocessing functions should all take the same input arguments,
% currently: unpack_dir, demodir

% Initialise empty struct with correct fields
% This to enforce additions to have the correct fields as well
datasets = struct('ID', {}, 'URL', {}, 'filetype', {}, 'postprocessing', {});

datasets(end+1) = struct( ...
    'ID', 'aa_demo', ...
    'URL', 'https://files.osf.io/v1/resources/umhtq/providers/osfstorage/5b7465680b87a00018c6a76f', ...
    'filetype', '.tar.gz', ...
    'postprocessing', @post_aa_demo);
datasets(end+1) = struct( ...
    'ID', 'ds000114', ...
    'URL', 'https://files.osf.io/v1/resources/9q7dv/providers/osfstorage/57e549f9b83f6901d457d162', ...
    'filetype', '.tar', ...
    'postprocessing', @post_ds000114);

end

function post_aa_demo(unpack_dir, demodir)
    movefile(fullfile(unpack_dir, 'aa_demo', '*'), demodir);
end

function post_ds000114(unpack_dir, demodir)
    movefile(fullfile(unpack_dir, 'ds114_test2', '*'), demodir);
end