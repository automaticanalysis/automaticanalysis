function aa_close(varargin)

switch nargin
    case 0
        aap.options.restorepath = 0;
        tasks = {'restorepath' 'restorewarnings' 'killjobs' 'clear'};
    case 1
        aap = varargin{1};
        tasks = {'restorepath' 'restorewarnings' 'killjobs' 'clear'};
    otherwise
        aap = varargin{1};
        tasks = varargin(2:end);
end

global aacache;

if cell_index(tasks,'restorepath') && ...
        (~isfield(aap.options,'restorepath') ||... % not specified
        aap.options.restorepath) % specified and enabled
    if isstruct(aacache) && isfield(aacache,'path')
        if isfield(aacache.path,'bcp_path'), path(aacache.path.bcp_path); end
        if isfield(aacache.path,'bcp_shellpath'), setenv('PATH', aacache.path.bcp_shellpath); end
    end
end

% restore warnings
if cell_index(tasks,'restorewarnings') && isstruct(aacache) && isfield(aacache,'warnings')
    for w = aacache.warnings
        warning(w);
    end
end

% kill running jobs
if cell_index(tasks,'killjobs')
    global taskqueue
    if isobject(taskqueue), taskqueue.close; end
end

% cleanup
if cell_index(tasks,'clear')
    clear global aa aacache aaparallel aaworker defaults taskqueue localtaskqueue;
end
end

