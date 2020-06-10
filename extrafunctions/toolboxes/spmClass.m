classdef spmClass < toolboxClass
    methods
        function obj = spmClass(path,varargin)
            defaultAddToPath = false;
            defaultKeepInPath = false;
            
            argParse = inputParser;
            argParse.addRequired('path',@ischar);
            argParse.addParameter('doAddToPath',defaultAddToPath,@(x) islogical(x) || isnumeric(x));
            argParse.addParameter('doKeepInPath',defaultKeepInPath,@(x) islogical(x) || isnumeric(x));
            argParse.parse(path,varargin{:});
            
            obj = obj@toolboxClass(argParse.Results.path,argParse.Results.doAddToPath,argParse.Results.doKeepInPath,{});
        end
        
        function load(obj)
            addpath(obj.toolPath);
            spm_jobman('initcfg');
            
            load@toolboxClass(obj)
        end
        
    end
end