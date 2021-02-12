function s_out = struct_update(varargin)

argParse = inputParser;
argParse.addRequired('s_in',@isstruct);
argParse.addRequired('s_upd',@isstruct);
argParse.addParameter('UpdateOnly',true,@(x)islogical(x)||isnumeric(x));
argParse.parse(varargin{:});
s_in = argParse.Results.s_in;
s_upd = argParse.Results.s_upd;

% Remove outdated fields
s_out = rmfield(s_in, intersect(fieldnames(s_in), fieldnames(s_upd)));
if argParse.Results.UpdateOnly, s_upd = rmfield(s_upd, setdiff(fieldnames(s_upd), fieldnames(s_in))); end

% Merge structs
s_out = cell2struct(...
    [struct2cell(s_out); struct2cell(s_upd)],...
    [fieldnames(s_out); fieldnames(s_upd)],...
    1);