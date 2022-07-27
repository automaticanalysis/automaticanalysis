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
            defaultRequiredPlugins = {};
            defaultKeepGUI = false;
            
            argParse = inputParser;
            argParse.addRequired('path',@ischar);
            argParse.addParameter('name','',@ischar);
            argParse.addParameter('doAddToPath',defaultAddToPath,@(x) islogical(x) || isnumeric(x));
            argParse.addParameter('requiredPlugins',defaultRequiredPlugins,@(x) iscellstr(x) || isstruct(x));
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
            
            obj = obj@toolboxClass(argParse.Results.name,argParse.Results.path,argParse.Results.doAddToPath,vars);
            
            obj.requiredPlugins = argParse.Results.requiredPlugins;
        end
        
        function load(obj,keepWorkspace)
            % backward compatibility
            if ~isstruct(obj.requiredPlugins)
                obj.requiredPlugins = cell2struct([obj.requiredPlugins;repmat({[]},2,10)],{'name','version','doPostprocess'});
            end
            
            if nargin < 2, keepWorkspace = false; end
            addpath(obj.toolPath);
            is_new_plugin = false;
            eeglab;
            if ~obj.showGUI, set(gcf,'visible','off'); end
            obj.plugins = evalin('base','PLUGINLIST');
            
            plPost = {};
            pllist = plugin_getweb('', obj.plugins, 'newlist');
            for p = reshape(pllist,1,[])
                if any(strcmp({obj.requiredPlugins.name}, p.name)) && ~p.installed
                    plInf = obj.requiredPlugins(strcmp({obj.requiredPlugins.name},p.name));
                    if ~isempty(plInf.version)
                        if isnumeric(plInf.version), plInf.version = num2str(plInf.version); end
                        p.version = plInf.version;
                        p.zip = spm_file(p.zip,'basename',[p.name p.version]);
                        p.size = 1; % force install
                    end
                    if ~isempty(plInf.doPostprocess) && plInf.doPostprocess, plPost(end+1) = {p.name}; end
                    
                    plugin_install(p.zip, p.name, p.version, p.size);
                    
                    is_new_plugin = true;
                end
            end
            if is_new_plugin, obj.load; end

            % postprocess 
            for p = reshape(plPost,1,[])
                aas_log([],false,['INFO: post-process installtion of ' p{1}]);
                obj.(['postprocess_' p{1}])();
            end
            
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
    
    methods (Access = private)
        function postprocess_AMICA(obj)
            [~,~,pl] = plugin_status('AMICA');
            plPath = fullfile(obj.toolPath,'plugins',pl.foldername);
            if isunix
               fileattrib(fullfile(plPath,'amica15ex'),'+x'); 
               fileattrib(fullfile(plPath,'amica15ub'),'+x'); 
               % shadow amica15ex with amica15ub (N.B.: revert if you work on the Expanse Supercomputer)
               movefile(fullfile(plPath,'amica15ex'),fullfile(plPath,'bcp_amica15ex'));
               system(sprintf('ln -s %s %s',fullfile(plPath,'amica15ub'),fullfile(plPath,'amica15ex')));
            end
        end
    end
end