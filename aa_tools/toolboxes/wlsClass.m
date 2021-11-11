% aa class interface for robust weighted least squares toolbox
%
% see:
% Diedrichsen, J., & Shadmehr, R. (2005). Detecting and adjusting for artifacts in fMRI time series data. Neuroimage, 27(3), 624-634
%

classdef wlsClass < toolboxClass
    
    properties (Access = protected)
        hGUI = []% GUI handles
    end
    
    methods
        
        function obj = wlsClass(path,varargin)
            
            defaultAddToPath = false;
            defaultKeepInPath = false;
            
            argParse = inputParser;
            argParse.addRequired('path', @ischar);
            argParse.addParameter('name','', @ischar);
            argParse.addParameter('doAddToPath', defaultAddToPath,@(x) islogical(x) || isnumeric(x));
            argParse.parse(path,varargin{:});
            
            obj = obj@toolboxClass(argParse.Results.name,argParse.Results.path,argParse.Results.doAddToPath,{});
        end
        
        function load(obj)
            addpath(obj.toolPath);
            load@toolboxClass(obj)
        end
        
    end
    
end