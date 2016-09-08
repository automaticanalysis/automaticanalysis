function aa_close(varargin)

if nargin
    tasks = varargin;
else
    tasks = {'restorepath' 'restorewarnings' 'killjobs' 'clear'};
end

% restore path
if cell_index(tasks,'restorepath') && ...
        (~isfield(aap.options,'restorepath') ||... % not specified
        aap.options.restorepath) % specified and enabled
    global aacache;
    if isstruct(aacache) && isfield(aacache,'bcp_path')
        path(aacache.bcp_path);
    end
    if isstruct(aacache) && isfield(aacache,'bcp_shellpath')
        setenv('PATH', aacache.bcp_shellpath);
    end
end

% restore warnings
if cell_index(tasks,'restorewarnings') && isfield(aacache,'warnings')
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

