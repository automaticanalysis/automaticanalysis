classdef hcpwbClass < toolboxClass
    properties
        binaryDir
        templateDir
    end
    
    properties (Access = protected)
        hGUI = []% GUI handles
    end
    
    methods
        function obj = hcpwbClass(path,varargin)
            defaultAddToPath = false;
            defaultKeepInPath = false;
            defaultTemplateDir = '';
            
            argParse = inputParser;
            argParse.addRequired('path',@ischar);
            argParse.addParameter('name','',@ischar);
            argParse.addParameter('doAddToPath',defaultAddToPath,@(x) islogical(x) || isnumeric(x));
            argParse.addParameter('doKeepInPath',defaultKeepInPath,@(x) islogical(x) || isnumeric(x));
            argParse.addParameter('templateDir',defaultTemplateDir,@ischar);
            argParse.parse(path,varargin{:});
            
            obj = obj@toolboxClass(argParse.Results.name,argParse.Results.path,argParse.Results.doAddToPath,argParse.Results.doKeepInPath,{});
            
            obj.templateDir = argParse.Results.templateDir;
        end
        
        function load(obj)
            % specify binary folder
            if isunix
                assert(~isempty(strfind(getenv('ARCH'),'64')),'64-bit OS required');
                if exist('/etc/os-release','file')
                    lines = regexp(fileread('/etc/os-release'), '\n', 'split');
                    lines(cellfun(@isempty,lines)) = [];
                    dat = regexp(lines, '\=', 'split');
                    ID=dat{cellfun(@(x) strcmp(x{1},'ID'),dat)}{2};
                    switch ID
                        case {'"centos"' '"rhel"' '"fedora"'}
                            obj.binaryDir = 'bin_rh_linux64';
                        otherwise
                            error('OS %s not supported', ID)
                    end
                else
                    error('OS not supported')
                end
            else
                error('/etc/os-release not found -> cannot determine *nix release')
            end
            
            % checkc templatedir
            if isempty(obj.templateDir) || ~exist(obj.templateDir,'dir') || numel(dir(obj.templateDir)) == 2
               warning('template directory is not found. You can obtain it from https://github.com/Washington-University/HCPpipelines (subfolder global/templates)');
            end
            
            load@toolboxClass(obj)
        end
    end
end