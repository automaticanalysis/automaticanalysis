function aaq_matlabpool(varargin)
switch varargin{1}
    case 'close'
        try matlabpool('close'); catch, delete(gcp('nocreate')); end
    otherwise
        try matlabpool(varargin{:}); catch, parpool(varargin{:}); end
end
end