classdef eeglabClass < toolboxClass
    properties
        req_plugins = {}
    end
    
    properties (Access = private)
        plugins = []
    end
    
    properties (Dependent)
        dipfit_path
    end
    
    methods
        function obj = eeglabClass(path,req_plugins)
            obj = obj@toolboxClass(path,true);
            if nargin > 1, obj.req_plugins = req_plugins; end
        end
        
        function init(obj)
            is_new_plugin = false;
            eeglab;
            close(gcf);
            obj.plugins = evalin('base','PLUGINLIST');
            evalin('base','clear ALLCOM ALLEEG CURRENTSET CURRENTSTUDY EEG eeglabUpdater globalvars LASTCOM PLUGINLIST STUDY')
            
            for p = plugin_getweb('', obj.plugins, 'newlist')
                if any(strcmp(obj.req_plugins, p.name)) && ~p.installed
                    plugin_install(p.zip, p.name, p.version, p.size);
                    is_new_plugin = true;
                end
            end
            if is_new_plugin, obj.init; end
            
            init@toolboxClass(obj)
        end
        
        function val = get.dipfit_path(obj)
            [~,~,pl] = plugin_status('dipfit');
            val = fullfile(obj.tool_path,'plugins',pl.foldername);
        end
    end
end