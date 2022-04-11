classdef spmClass < toolboxClass
    properties (Access = protected)
        hGUI = []% GUI handles
    end
    
    methods
        function obj = spmClass(path,varargin)
            defaultAddToPath = false;
            
            argParse = inputParser;
            argParse.addRequired('path',@ischar);
            argParse.addParameter('name','',@ischar);
            argParse.addParameter('doAddToPath',defaultAddToPath,@(x) islogical(x) || isnumeric(x));
            argParse.parse(path,varargin{:});
            
            obj = obj@toolboxClass(argParse.Results.name,argParse.Results.path,argParse.Results.doAddToPath,{});
            
            obj.addToolbox(fieldtripClass(fullfile(obj.toolPath,'external','fieldtrip'),'name','fieldtrip'));
            
            obj.collections(1).name = 'meeg';
            obj.collections(1).path = {...
                    'external/bemcp'...
                    'external/ctf'...
                    'external/eeprobe'...
                    'external/mne'...
                    'external/yokogawa_meg_reader'...
                    'toolbox/dcm_meeg'...
                    'toolbox/spectral'...
                    'toolbox/Neural_Models'...
                    'toolbox/MEEGtools'...
                    };
            obj.collections(1).toolbox = {'fieldtrip'};
        end
        
        function load(obj)
            addpath(obj.toolPath);
            
            load@toolboxClass(obj)
            
            spm_jobman('initcfg');
            
            % Update stored path
            p = split(path,pathsep);                
            obj.toolInPath = cat(1,...
                obj.toolInPath,...
                setdiff(p(cellfun(@(x) contains(x,obj.toolPath), p)),obj.toolInPath)...
                );
        end
    end
end