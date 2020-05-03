classdef eeglabClass < toolboxClass
    properties
        requiredPlugins = {}
    end
    
    properties (Access = private)
        plugins = []
    end
    
    properties (Dependent)
        dipfitPath
    end
    
    methods
        function obj = eeglabClass(path,varargin)
            defaultKeepInPath = false;
            defaultRequiredPlugins = {};
            
            argParse = inputParser;
            argParse.addRequired('path',@ischar);
            argParse.addParameter('requiredPlugins',defaultRequiredPlugins,@iscellstr);
            argParse.addParameter('doKeepInPath',defaultKeepInPath);
            argParse.parse(path,varargin{:});
            
            obj = obj@toolboxClass(argParse.Results.path,true,argParse.Results.doKeepInPath);
            
            obj.requiredPlugins = argParse.Results.requiredPlugins;
        end
        
        function load(obj)
            is_new_plugin = false;
            eeglab;
            close(gcf);
            obj.plugins = evalin('base','PLUGINLIST');
            evalin('base','clear ALLCOM ALLEEG CURRENTSET CURRENTSTUDY EEG eeglabUpdater globalvars LASTCOM PLUGINLIST STUDY')
            
            for p = plugin_getweb('', obj.plugins, 'newlist')
                if any(strcmp(obj.requiredPlugins, p.name)) && ~p.installed
                    plugin_install(p.zip, p.name, p.version, p.size);
                    is_new_plugin = true;
                end
            end
            if is_new_plugin, obj.init; end
            
            load@toolboxClass(obj)
        end
        
        function val = get.dipfitPath(obj)
            [~,~,pl] = plugin_status('dipfit');
            val = fullfile(obj.tool_path,'plugins',pl.foldername);
        end
    end
end