classdef spmClass < toolboxClass
    properties (Access=private)
        externalSpec = [...
            struct('name','meeg',...
            'path','external/bemcp:external/ctf:external/eeprobe:external/mne:external/yokogawa_meg_reader:toolbox/dcm_meeg:toolbox/spectral:toolbox/Neural_Models:toolbox/MEEGtools'...
            ) ...
            ]
    end
    
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
        
        function addExternal(obj,toolbox)
            if obj.pStatus < obj.CONST_STATUS.loaded
                warning('SPM is not loaded')
                return
            end
            iTbx = strcmp(obj.externalSpec.name,toolbox);
            if ~any(iTbx)
                warning('external %s not specified',toolbox);
                return
            end
            for p = split(obj.externalSpec(iTbx).path,':')'
                currp = fullfile(obj.toolPath,strrep(p{1},'/',filesep));
                addpath(currp);
                obj.toolInPath = vertcat(obj.toolInPath,currp);
            end
        end
        
        function rmExternal(obj,toolbox)
            iTbx = strcmp(obj.externalSpec.name,toolbox);
            if ~any(iTbx)
                warning('external %s not specified',toolbox);
                return
            end
            for p = split(obj.externalSpec(iTbx).path,':')'
                currp = fullfile(obj.toolPath,strrep(p{1},'/',filesep));
                rmpath(currp);
                obj.toolInPath(strcmp(obj.toolInPath,currp)) = [];
            end
        end
    end
end