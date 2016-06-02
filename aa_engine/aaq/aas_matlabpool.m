function out = aas_matlabpool(varargin)
switch varargin{1}
    case 'getcurrent'
        try out = gcp('nocreate'); catch
            warning('Not supported by MATLAB %s',sprintf('%d.',getmatlabversion));
            out = [];
        end
    case 'isopen'
        try out = ~isempty(gcp('nocreate')); catch, out = matlabpool('size') > 0; end
    case 'close'
        try 
            p = gcp('nocreate');
            if ~isempty(p), delete(p); end
        catch
            try
                matlabpool('close'); 
            catch E
                if strcmp(E.identifier, 'parallel:convenience:InteractiveSessionNotRunning')
                    warning('No open pool is found');
                else
                    rethrow(E);
                end
            end
        end
    otherwise % start
        try out = parpool(varargin{:}); catch
            matlabpool(varargin{:}); 
            out = [];
        end
end
end