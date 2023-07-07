% Automatic analysis - changes root directory
% This function changes the root directory by changing the root path in
% each aap structure within an analysis results directory.
% Steps:
% 1. Move your results (ie., folder with the name analysisid; aap.directory_conventions.analysisid) to the new directory. 
% 2. Call this function by giving the old root directory, the new root directory (where you have moved the files) and the analysis_id (name of the folder where your results are)
% Parameters:
%       old_root: corresponds to the old root directory (aap.acq_details.root)
%       new_root: The new root directory (the directory you moved your analysisid).
%       analysis_id: your aap.directory_conventions.analysisid
% Tamer Gezici, Bilkent 2023 -- initial implementation

function aa_change_root(old_root,new_root,analysis_id)
search_path = fullfile(new_root,analysis_id); % New root directory + results
[files,dirs] = spm_select('FPListRec', search_path, '.*aap_parameters.*\.mat$'); % Find aap parameters
files = cellstr(files);
for i=1:length(files)
    load(files{i},"aap");
    old_path = aap.acq_details.root;
    new_path = strrep(old_path,old_root,new_root);
    aap.acq_details.root = new_path;
    aap.aap_beforeuserchanges.acq_details.root = new_path;
    save(files{i},"aap");
end
end

