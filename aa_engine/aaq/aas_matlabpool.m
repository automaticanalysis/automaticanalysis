function out = aas_matlabpool(varargin)
% AAS_MATLABPOOL performs actions on the current matlab pool:
%   out = aas_matlabpool('getcurrent') returns the current pool object,
%       which may be empty (no pool will be started)
%   out = aas_matlabpool('isopen') returns a logical indicating whether a 
%       pool is open 
%   out = aas_matlabpool('close') closes the current pool if it exists and
%       returns an empty pool object
%   out = aas_matlabpool(varargin) where varargin is none of the options 
%       above, calls parpool; varargin in this case must be a resource, or 
%       a poolsize, or both, plus name-value pairs (see parpool). The pool
%       object will be returned.
% 
% If the request fails, [] will be returned.

% default value upon failure
out = [];
try
    switch varargin{1}
        case 'getcurrent'
            out = gcp('nocreate');
        case 'isopen'
            out = ~isempty(gcp('nocreate'));
        case 'close'
            delete(gcp('nocreate'));
            % return empty pool project rather than []
            out = gcp('nocreate');
        otherwise % start
            out = parpool(varargin{:});
    end
catch ME
    warning(ME.identifier, 'An error was caught handling the parallel pool:\n %s ', ME.message);
end

end