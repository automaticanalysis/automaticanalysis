function is_absolute = is_absolute_path(fpath)
% is_absolute = is_absolute_path(fpath)
% Return boolean whether input fpath is an absolute path

% make sure fileseperators are right
fpath = fullfile(fpath, '');

% check if path is absolute
is_absolute = true;
if ispc()
    if length(fpath) < 2 || (fpath(2) ~= ':' && ~strcmpi(repmat(filesep,1,2), fpath(1:2)))
        is_absolute = false;
    end
elseif isunix()
    if ~(fpath(1) == filesep || fpath(1) == '~')
        is_absolute = false;
    end
else
    error('Unsupported operating system');
end

end