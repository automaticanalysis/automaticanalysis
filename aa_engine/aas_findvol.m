function subjdir = aas_findvol(aap,subj)

% convert to new format
if ~isstruct(aap.directory_conventions.rawdatadir)
    SEARCHPATH{1} = aap.directory_conventions.rawdatadir;
elseif ~iscell(aap.directory_conventions.rawdatadir.paths)
    SEARCHPATH{1} = aap.directory_conventions.rawdatadir.paths;
else
    SEARCHPATH = aap.directory_conventions.rawdatadir.paths;
end

isFound = false;
for i = 1:numel(SEARCHPATH)
	if isnumeric(subj) % search among subjects already added
		subjdir=fullfile(SEARCHPATH{i},aap.acq_details.subjects(subj).mriname);
	elseif ischar(subj) % custom search
		subjdir=fullfile(SEARCHPATH{i},subj);
    else
        aas_log(aap,1,sprintf('Input must be either integer or string!'));
	end
    if exist(subjdir,'dir')
        isFound = true;
        break; 
    end
end

if ~isFound, subjdir = ''; end
    