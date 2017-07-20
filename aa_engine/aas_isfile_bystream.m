function isanyfile = aas_isfile_bystream(aap,varargin)
aap.options.verbose = -1;
isanyfile = ~isempty(aas_getfiles_bystream(aap,varargin{:}));