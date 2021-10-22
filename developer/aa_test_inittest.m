function aap = aa_test_inittest(testpath, parameterfile, deleteprevious, wheretoprocess)

% this is a convenience function called by testing scripts to setup
% data and results directories and required processing options

if isempty(parameterfile)
    aap = aarecipe([testpath '.xml']);
else
    aap = aarecipe(parameterfile,[testpath '.xml']);
end

temp = strsplit(testpath,filesep); temp = strsplit(temp{end},'_');

aap.directory_conventions.analysisid = [ temp{2} '_' temp{3} ];
aap.directory_conventions.rawdatadir = fullfile(aap.directory_conventions.rawdatadir,temp{2});

anadir = fullfile(aap.acq_details.root, aap.directory_conventions.analysisid);
fprintf('Saving results in: %s\n', anadir);
if exist(anadir,'dir') && deleteprevious
    fprintf('Removing previous results...');
    rmdir(anadir,'s');
    fprintf('Done\n');    
end

aap.options.wheretoprocess = wheretoprocess;