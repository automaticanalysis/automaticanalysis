function aap = aa_test_inittest(testpath,deleteprevious)

aap = aarecipe([testpath '.xml']);

% -------------------------------------------------------------------------
% results and data directory specification
% -------------------------------------------------------------------------

temp = strsplit(testpath,filesep); temp = strsplit(temp{end},'_');
aap.directory_conventions.analysisid = [ temp{2} '_' temp{3} ];

anadir = fullfile(aap.acq_details.root, aap.directory_conventions.analysisid);
fprintf('Saving results in: %s\n', anadir);
if exist(anadir,'dir') && deleteprevious
    fprintf('Removing previous results...');
    rmdir(anadir,'s');
    fprintf('Done\n');    
end

if contains(aap.directory_conventions.rawdatadir,':') % more then one directories
    demodir = regexp(aap.directory_conventions.rawdatadir,'(?<=:)[a-zA-Z0-9_\/]*aa_demo(?:)','match');
else % one directory
    demodir = regexp(aap.directory_conventions.rawdatadir,'[a-zA-Z0-9_\/]*aa_demo(?:)','match');
end
if isempty(demodir), aas_log([],true,'aap.directory_conventions.rawdatadir contains no aa_demo directory'); end   
aap.directory_conventions.rawdatadir = fullfile(demodir{1},temp{2});