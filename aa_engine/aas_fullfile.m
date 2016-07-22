% Like full file, except that input filepaths is a cell array. basepath is
% prepended (using fullfile) to every item in filepaths. Result is a cell
% array

function fullfilepaths=aas_fullfile(basepath,filepaths)

for fileind=1:length(filepaths)
    fullfilepaths{fileind}=fullfile(basepath,filepaths{fileind});
end;
