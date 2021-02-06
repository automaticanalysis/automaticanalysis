classdef toolboxClass < handle
    properties
        name = ''
        toolPath = ''
        autoLoad = false % (re-)load sub-toolbox with its parent
        keepInPath = false % keep the toolbox in the path after deletion
        showGUI = false % show GUI upon load
    end
    
    properties (Dependent)
        status
    end
    
    properties (Access = protected)
        CONST_STATUS = struct(...
            'undefined', -1,...
            'defined', 0,...
            'unloaded', 1,...
            'loaded', 2 ...            
        )
        pStatus = -1
        
        toolInPath = {}
        
        workspace = struct
        
        collections = struct('name',{},'path',{},'toolbox',{})
        
        toolboxes = cell(1,0)
    end
    
    properties (Access = protected, Abstract)
        hGUI % GUI handles
    end
    
    methods
        function obj = toolboxClass(name,path,doAddToPath,doKeepInPath,workspaceVariableNames)
            obj.name = name;
            obj.toolPath = path;
            if doAddToPath
                addpath(obj.toolPath); 
                obj.toolInPath = cellstr(obj.toolPath);
            end
            obj.keepInPath = doKeepInPath;
            obj.pStatus = obj.CONST_STATUS.defined;
            
            for v = workspaceVariableNames
                obj.workspace.(v{1}) = [];
            end
        end
        
        function val = get.status(obj)
            for f = fieldnames(obj.CONST_STATUS)'
                if obj.CONST_STATUS.(f{1}) == obj.pStatus
                    val = f{1};
                    break;
                end
            end
        end
        
        function load(obj,keepWorkspace)
            if nargin < 2, keepWorkspace = false; end
            if obj.pStatus < obj.CONST_STATUS.loaded
                p = split(path,pathsep);                
                obj.toolInPath = p(cellfun(@(x) ~isempty(strfind(x,obj.toolPath)), p));
                obj.pStatus = obj.CONST_STATUS.loaded;
            end
            for v = fieldnames(obj.workspace)'
                if ~any(strcmp(evalin('base','who'),v{1}))
                    obj.workspace = rmfield(obj.workspace,v{1});
                    continue;
                end
                obj.workspace.(v{1}) = evalin('base',v{1});
                if ~keepWorkspace, evalin('base',['clear ' v{1}]); end
            end
        end
        
        function close(obj)
            % remove from path
            if ~obj.keepInPath, obj.unload; end
            
            % close sub-toolboxes
            for t = obj.toolboxes
                t{1}.close;
            end
            
            % GUI and workspace
            for h = obj.hGUI, close(h); end
            for v = fieldnames(obj.workspace)'
                obj.workspace.(v{1}) = [];
            end
            
            obj.pStatus = obj.CONST_STATUS.undefined;
        end
        
        function reload(obj,loadWorkspace)
            % default
            if nargin < 2, loadWorkspace = false; end
            
            % re-add to path
            if obj.pStatus < obj.CONST_STATUS.loaded
                addpath(sprintf(['%s' pathsep],obj.toolInPath{:}))
                obj.pStatus = obj.CONST_STATUS.loaded;
            end
            
            % reload sub-toolboxes
            for t = obj.toolboxes
                if t{1}.autoLoad, t{1}.reload; end
            end
            
            % GUI and workspace
            if obj.showGUI, for h = obj.hGUI, set(h,'visible','on'); end; end
            if loadWorkspace
                for v = fieldnames(obj.workspace)'
                    assignin('base', v{1}, obj.workspace.(v{1}));
                end
            end
        end
        
        function unload(obj,updateWorkspace)
            % default
            if nargin < 2, updateWorkspace = false; end
            
            % remove from path
            if obj.pStatus > obj.CONST_STATUS.unloaded
                warning('%s''s folders (and subfolders) will be removed from the MATLAB path',class(obj));
                rmpath(sprintf(['%s' pathsep],obj.toolInPath{:}))
                obj.pStatus = obj.CONST_STATUS.unloaded;
            end
            
            % unload sub-toolboxes
            for t = obj.toolboxes
                t{1}.unload;
            end
            
            % GUI and workspace
            if obj.showGUI, for h = obj.hGUI, set(h,'visible','off'); end; end
            for v = fieldnames(obj.workspace)'
                if updateWorkspace, obj.workspace.(v{1}) = evalin('base',v{1}); end
                evalin('base',['clear ' v{1}]);
            end
        end
        
        %% sub-toolboxes
        function addToolbox(obj,tb)
            if isempty(tb.name)
                warning('sub-toolbox MUST have a name');
                return
            end
            isTb = cellfun(@(t) strcmp(tb.name,t.name), obj.toolboxes);
            if isTb
                warning('toolbox %s is already added as a sub-toolbox', tb.name);
            else
                if tb.autoLoad, tb.load; end
                obj.toolboxes{end+1} = tb;                
            end
        end
        
        function doToolbox(obj,tbname,task)
            if ~any(strcmp(methods('toolboxClass'),task))
                warning('unsupported task: %s',task)
                return
            end
            
            isTb = cellfun(@(t) strcmp(tbname,t.name), obj.toolboxes);
            if ~isTb
                warning('toolbox %s is not a sub-toolbox', tbname);
            else
                obj.toolboxes{isTb}.(task);
            end
        end
        
        function rmToolbox(obj,tbname)
            isTb = cellfun(@(t) strcmp(tbname,t.name), obj.toolboxes);
            if ~isTb
                warning('toolbox %s is not a sub-toolbox', tbname);
            else
                tb = obj.toolboxes{isTb};
                tb.close;
                obj.toolboxes(isTb) = [];
            end
        end
        
        function setAutoLoad(obj)
            obj.autoLoad = true;
        end
        function unsetAutoLoad(obj)
            obj.autoLoad = false;
        end
        
        %% collections
         function addCollection(obj,collection)
            if obj.pStatus < obj.CONST_STATUS.loaded
                warning('toolbox is not loaded')
                return
            end
            iC = strcmp(obj.collections.name,collection);
            if ~any(iC)
                warning('external %s not specified',collection);
                return
            end
            for p = obj.collections(iC).path
                currp = fullfile(obj.toolPath,strrep(p{1},'/',filesep));
                addpath(currp);
                obj.toolInPath = vertcat(obj.toolInPath,currp);
            end
            for t = obj.collections(iC).toolbox
                obj.doToolbox(t{1},'load');
            end
        end
        
        function rmCollection(obj,collection)
            iC = strcmp(obj.collections.name,collection);
            if ~any(iC)
                warning('external %s not specified',collection);
                return
            end
            for p = obj.collections(iC).path
                currp = fullfile(obj.toolPath,strrep(p{1},'/',filesep));
                rmpath(currp);
                obj.toolInPath(strcmp(obj.toolInPath,currp)) = [];
            end
            for t = obj.collections(iC).toolbox
                obj.doToolbox(t{1},'unload');
            end
        end
    end
end