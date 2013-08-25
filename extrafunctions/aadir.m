% return the path to current aa4 directory
% d = aadir;
function d = aadir;

% start aa if it isn't already running
if ~exist('aas_addsubject') == 2
    aa_ver4_devel;
end

% find the parent dir
d = fileparts(fileparts(which('aas_addsubject')));
