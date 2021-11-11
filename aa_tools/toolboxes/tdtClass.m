classdef tdtClass < toolboxClass
    properties
        defaults = []% defaults
    end
    
    properties (Access = protected)
        hGUI = []% GUI handles
    end
    
    methods
        function obj = tdtClass(path,varargin)
            defaultAddToPath = false;
            
            argParse = inputParser;
            argParse.addRequired('path',@ischar);
            argParse.addParameter('name','',@ischar);
            argParse.addParameter('doAddToPath',defaultAddToPath,@(x) islogical(x) || isnumeric(x));
            argParse.parse(path,varargin{:});
            
            obj = obj@toolboxClass(argParse.Results.name,argParse.Results.path,argParse.Results.doAddToPath,{});
        end
        
        function load(obj)
            addpath(obj.toolPath);
            obj.defaults = decoding_defaults;
            
            load@toolboxClass(obj)
        end
    end
end