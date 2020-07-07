classdef fieldtripClass < toolboxClass
    methods
        function obj = fieldtripClass(path,varargin)
            defaultAddToPath = false;
            defaultKeepInPath = false;
            
            argParse = inputParser;
            argParse.addRequired('path',@ischar);
            argParse.addParameter('name','',@ischar);
            argParse.addParameter('doAddToPath',defaultAddToPath,@(x) islogical(x) || isnumeric(x));
            argParse.addParameter('doKeepInPath',defaultKeepInPath,@(x) islogical(x) || isnumeric(x));
            argParse.parse(path,varargin{:});
            
            obj = obj@toolboxClass(argParse.Results.name,argParse.Results.path,argParse.Results.doAddToPath,argParse.Results.doKeepInPath,{});
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
        
        function addExternal(obj,toolbox)
            if obj.pStatus < obj.CONST_STATUS.loaded
                warning('Fieldtrip is not loaded')
                return
            end
            tbpath = '';
            lastwarn('')
            ft_hastoolbox(toolbox,1);
            if ~isempty(strfind(lastwarn,toolbox))
                tbpath = strsplit(lastwarn); tbpath = tbpath{2};
            end
            p = split(path,pathsep);                
            obj.toolInPath = vertcat(obj.toolInPath, p(cellfun(@(x) ~isempty(strfind(x,tbpath)), p)));
            
            if ~isempty(strfind(toolbox,'spm'))
                spmVer = spm('ver'); 
                global ft_default
                ft_default.trackcallinfo = 'no';
                ft_default.showcallinfo = 'no';
                if ~isempty(spmVer), ft_default.spmversion = lower(spmVer); end
                fprintf('INFO: SPM version %s is set\n',spmVer);
            end
        end
        
        function rmExternal(obj,toolbox)
            indtbp = cellfun(@(x) ~isempty(strfind(x,toolbox)), obj.toolInPath);
            p = obj.toolInPath(indtbp);
            rmpath(sprintf(['%s' pathsep],p{:}))
            obj.toolInPath(indtbp) = [];
        end
    end
end