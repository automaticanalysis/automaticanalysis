
classdef eeglabClass < toolboxClass
    properties
        requiredPlugins = {}
    end
    
    properties (Access = private)
        plugins = []
    end
    
    properties (Access = protected, Dependent)
        hGUI% GUI handles
    end
    
    properties (Dependent)
        dipfitPath
    end
    
    methods
        function obj = eeglabClass(path,varargin)
            defaultAddToPath = true;
            defaultKeepInPath = false;
            defaultRequiredPlugins = {};
            defaultKeepGUI = false;
            
            argParse = inputParser;
            argParse.addRequired('path',@ischar);
            argParse.addParameter('name','',@ischar);
            argParse.addParameter('doAddToPath',defaultAddToPath,@(x) islogical(x) || isnumeric(x));
            argParse.addParameter('doKeepInPath',defaultKeepInPath,@(x) islogical(x) || isnumeric(x));
            argParse.addParameter('requiredPlugins',defaultRequiredPlugins,@iscellstr);
            argParse.addParameter('doKeepGUI',defaultKeepGUI,@(x) islogical(x) || isnumeric(x));
            argParse.parse(path,varargin{:});
            
            vars = {...
                '{"name": "ALLCOM"}'...
                '{"name": "ALLEEG"}'...
                '{"name": "CURRENTSET"}'...
                '{"name": "CURRENTSTUDY"}'...
                '{"name": "eeglabUpdater"}'...
                '{"name": "globalvars"}'...
                '{"name": "LASTCOM"}'...
                '{"name": "PLUGINLIST"}'...
                '{"name": "STUDY"}'...
                };
            
            obj = obj@toolboxClass(argParse.Results.name,argParse.Results.path,argParse.Results.doAddToPath,argParse.Results.doKeepInPath,vars);
            
            obj.requiredPlugins = argParse.Results.requiredPlugins;
        end
        
        function load(obj,keepWorkspace)
            if nargin < 2, keepWorkspace = false; end
            addpath(obj.toolPath);
            is_new_plugin = false;
            eeglab;
            if ~obj.showGUI, set(gcf,'visible','off'); end
            obj.plugins = evalin('base','PLUGINLIST');
            
            pllist = plugin_getweb('', obj.plugins, 'newlist');
            for p = reshape(pllist,1,[])
                if any(strcmp(obj.requiredPlugins, p.name)) && ~p.installed
                    plugin_install(p.zip, p.name, p.version, p.size);
                    is_new_plugin = true;
                end
            end
            if is_new_plugin, obj.load; end
            
            load@toolboxClass(obj,keepWorkspace)
        end
        
        function val = get.hGUI(obj)
            if ~obj.showGUI, val = [];
            else
                val = findall(allchild(0),'Flat','Tag','EEGLAB');
            end
        end
        
        function val = get.dipfitPath(obj)
            [~,~,pl] = plugin_status('dipfit');
            val = fullfile(obj.toolPath,'plugins',pl.foldername);
        end
    end
end