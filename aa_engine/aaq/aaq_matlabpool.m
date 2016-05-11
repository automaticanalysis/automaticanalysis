function out = aaq_matlabpool(varargin)
switch varargin{1}
    case 'isopen'
        try out = matlabpool('size') > 0; catch, out = ~isempty(gcp('nocreate')); end
    case 'close'
        try matlabpool('close'); catch, delete(gcp('nocreate')); end
    otherwise
        try matlabpool(varargin{:}); catch, parpool(varargin{:}); end
end
end