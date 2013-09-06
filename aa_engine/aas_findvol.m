function subjdir = aas_findvol(aap,subj)

% convert to new format
if ~isstruct(aap.directory_conventions.rawdatadir)
    warning('off','MATLAB:warn_r14_stucture_assignment')
    aap.directory_conventions.rawdatadir.paths{1} = aap.directory_conventions.rawdatadir;
    warning('on','MATLAB:warn_r14_stucture_assignment')
end
SEARCHPATH = aap.directory_conventions.rawdatadir.paths;

isFound = false;
for i = 1:numel(SEARCHPATH)
    subjdir=fullfile(SEARCHPATH{i},aap.acq_details.subjects(subj).mriname);
    if exist(subjdir,'dir')
        isFound = true;
        break; 
    end
end

if ~isFound, subjdir = ''; end
    