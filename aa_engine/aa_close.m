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

if cell_index(tasks,'restorepath') && ...
        (~isstruct(aap) || ~isfield(aap.options,'restorepath') ||... % not specified
        aap.options.restorepath) % specified and enabled
    [s,p] = aas_cache_get(aap,'bcp_path','system'); if s, path(p); end
    [s,p] = aas_cache_get(aap,'bcp_shellpath','system'); if s, setenv('PATH', p); end
end

% restore warnings
if cell_index(tasks,'restorewarnings')
    [s,ws] = aas_cache_get(aap,'warnings','system');
    if s
       for w = ws
           warning(w);
       end
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

% close any SPM windows that may have been left up
close all

end

