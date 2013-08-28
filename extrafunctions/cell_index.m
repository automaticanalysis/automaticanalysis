% Locate string "str" within a cell of strings "lines"
% Returns 0 if none
% Optional parameters:
% 	nl:	starting index (default = 1)
% 	ni:	returns location only if "str" is a substring starting from "ni"th character (default = none --> any match)
% Tibor Auer MRC CBU Cambridge 2012-2013

function r = cell_index(varargin)
i = sub(varargin{:});
r = [];
while true
    r(end+1) = i;
    if i == numel(varargin{1}), break, end
    i = sub(varargin{:},i+1);    
    if ~i, break, end
end
end
    
function i = sub(varargin)
lines = varargin{1};
str = varargin{2};
if (nargin < 3) || isempty(varargin{3})
    nl = 1;
else
    nl = varargin{3};
end

for i = nl:numel(lines)
    sample = lines{i};
    if iscell(sample) || (size(sample,1) > 1)
        sample = [];
        continue; 
    end
    if size(sample,1) > size(sample,2), sample = sample';  end
    if ~ischar(sample), sample = num2str(sample); end

    sample = strfind(sample,str);
    if ~isempty(sample) && ((nargin < 4) || ((nargin > 3) && any(sample == varargin{4}))), break, end
end
if isempty(sample)
    i = 0;
end
end