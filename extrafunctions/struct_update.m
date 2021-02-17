function s_out = struct_update(varargin)

argParse = inputParser;
argParse.addRequired('s_in',@isstruct);
argParse.addRequired('s_upd',@isstruct);
argParse.addParameter('Mode','update',@ischar);
argParse.parse(varargin{:});
s_in = argParse.Results.s_in;
s_upd = argParse.Results.s_upd;

switch argParse.Results.Mode
    case 'update'
        % remove common fields from input
        s_out = rmfield(s_in, intersect(fieldnames(s_in), fieldnames(s_upd)));
        % remove missing fields from update
        s_upd = rmfield(s_upd, setdiff(fieldnames(s_upd), fieldnames(s_in)));
    case 'extend'
        s_out = s_in;
        % remove common fields from update
        s_upd = rmfield(s_upd, intersect(fieldnames(s_upd), fieldnames(s_in)));
    otherwise % update and extend
        % remove common fields from input
        s_out = rmfield(s_in, intersect(fieldnames(s_in), fieldnames(s_upd)));
end

% Merge structs
s_out = cell2struct(...
    [struct2cell(s_out); struct2cell(s_upd)],...
    [fieldnames(s_out); fieldnames(s_upd)],...
    1);