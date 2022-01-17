function [img, md5, streamdesc, domain] = aas_getfiles_bystream_multilevel(aap,varargin)
% Multilevel streamlocator: session/specified --> subject --> study --> lower ...
% Inputs:
%   aap
%   domain
%   index/indices
%   streamname
%   source: 'input'/'output' (optional, default: 'input')

source = 'input';
if strcmp(varargin{end},'input') || strcmp(varargin{end},'output')
    source = varargin{end}; 
    varargin(end) = [];
end
domain = varargin{1};
index = [varargin{2:end-1}];
streamname = varargin{end};

img = ''; d = 0;
verbose0 = aap.options.verbose; aap.options.verbose = -1; % ignore error

% parents
while isempty(img) && ~isempty(domain)
    d = d + 1;
    if d > 1, domain = aas_finddomain(aap,domain,'parent'); end
    if ~isempty(domain), [img, md5, streamdesc] = aas_getfiles_bystream(aap,domain,index(1:end-d+1),streamname,source); end
end

if isempty(img)
    % children
    domain = aas_finddomain(aap,varargin{1},'children');
    index = [index 1];
    while ~isempty(domain)
        img = cellfun(@(d) aas_getfiles_bystream(aap,d,index,streamname,source),domain,'UniformOutput',false);
        idx = find(cellfun(@(x) ~isempty(x), img),1,'first');
        if ~isempty(idx), img = img{idx}; break; end
        domain = aas_finddomain(aap,domain,'children');
        index = [index 1];
    end
end

aap.options.verbose = verbose0;

if isempty(img)
    aas_log(aap,true,sprintf('%s stream %s not found',source,streamname));
end
end

function outdomain = aas_finddomain(aap,indomain,direction,start)
% find children or parent domain depending on the direction
if nargin < 4 || isempty(start), start = {'directory_conventions' 'parallel_dependencies'}; end
outdomain = {};
if ~isempty(getfield(aap,start{:}))
    tree = getfield(aap,start{:});
    fn = fieldnames(tree);
    start0 = start;
    for fnind = 1:length(fn)
        if strcmp(fn{fnind},indomain)
            switch direction
                case 'parent'
                    if ~strcmp(start0{end},'parallel_dependencies')
                        outdomain = start0{end};
                    end
                case 'children'
                    if isstruct(tree.(char(indomain)))
                        outdomain = fieldnames(tree.(indomain));
                    end
            end
            return;
        else
            start = [start0 fn(fnind)];
            outdomain = aas_finddomain(aap,indomain,direction,start);
            if ~isempty(outdomain)
                return;
            end
        end
    end
end
end

function [label index] = aa_unique(vals)
try
    [label index] = unique(vals, 'legacy');
catch
    [label index] = unique(vals);
end
end