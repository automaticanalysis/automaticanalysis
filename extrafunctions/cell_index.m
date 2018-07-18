% Locate string "str" within a cell of strings "lines"
% Returns 0 if none
% Optional parameters:
% 	nl:	starting index (default = 1)
% 	ni:	returns location only if "str" is a substring starting from "ni"th character (default = none --> any match)
% Tibor Auer MRC CBU Cambridge 2012-2013

function r = cell_index(varargin)
nl = 1;
varargin{1}(cellfun(@(x) ~ischar(x), varargin{1})) = {''}; % ingore nonchar elements

if (nargin >= 3) && ~isempty(varargin{3}), nl = varargin{3}; end
if (nargin >= 4)
    r = find(cellfun(@(x) ~isempty(x) && (x == varargin{4}),strfind(varargin{1}(nl:end),varargin{2})));
else
    r = find(cellfun(@(x) ~isempty(x),strfind(varargin{1}(nl:end),varargin{2})));
end
if isempty(r), r = 0; end
r = r';
end