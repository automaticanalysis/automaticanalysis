classdef fieldtripClass < toolboxClass
    methods
        function obj = fieldtripClass(path,varargin)
            defaultAddToPath = false;
            defaultKeepInPath = false;
            
            argParse = inputParser;
            argParse.addRequired('path',@ischar);
            argParse.addParameter('doAddToPath',defaultAddToPath,@(x) islogical(x) || isnumeric(x));
            argParse.addParameter('doKeepInPath',defaultKeepInPath,@(x) islogical(x) || isnumeric(x));
            argParse.parse(path,varargin{:});
            
            obj = obj@toolboxClass(argParse.Results.path,argParse.Results.doAddToPath,argParse.Results.doKeepInPath);
        end
        
        function load(obj)
            addpath(obj.toolPath);
            ft_defaults
            spmVer = '';
            try spmVer = spm('ver'); catch, warning('SPM is not detected'); end
            if ~isempty(spmVer)
                global ft_default
                ft_default.trackcallinfo = 'no';
                ft_default.showcallinfo = 'no';
                if ~isempty(spmVer), ft_default.spmversion = lower(spmVer); end
            end
            
            load@toolboxClass(obj)
        end
        
        function close(obj)
            global ft_default
            ft_default = [];
            clear ft_default;
            close@toolboxClass(obj)
        end        
    end
end