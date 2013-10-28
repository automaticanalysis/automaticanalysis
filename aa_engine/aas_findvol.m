function subjdir = aas_findvol(aap,subj)
%% Changed to comma separated list format [RC]
if isstruct(aap.directory_conventions.rawdatadir)
    aas_log(aap,true,'Structure formate for rawdatadir no longer supported - use comma separated list');
end;

% Parse comma separated list
SEARCHPATH = textscan(aap.directory_conventions.rawdatadir,'%s','delimiter', ':');
SEARCHPATH = SEARCHPATH{1};

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
    