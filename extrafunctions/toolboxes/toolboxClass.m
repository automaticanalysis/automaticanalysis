classdef toolboxClass < handle
    properties
        toolPath = ''
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
        
        hGUI % GUI handles
    end
    
    methods
        function obj = toolboxClass(path,doAddToPath,doKeepInPath)
            obj.toolPath = path;
            if doAddToPath
                addpath(obj.toolPath); 
                obj.toolInPath = cellstr(obj.toolPath);
            end
            obj.keepInPath = doKeepInPath;
            obj.pStatus = obj.CONST_STATUS.defined;
        end
        
        function val = get.status(obj)
            for f = fieldnames(obj.CONST_STATUS)'
                if obj.CONST_STATUS.(f{1}) == obj.pStatus
                    val = f{1};
                    break;
                end
            end
        end
        
        function load(obj)
            if obj.pStatus < obj.CONST_STATUS.loaded
                p = split(path,pathsep);
                obj.toolInPath = p(cellfun(@(x) contains(x,obj.toolPath), p));
                obj.pStatus = obj.CONST_STATUS.loaded;
            end
        end
        
        function close(obj)
            if ~obj.keepInPath, obj.unload; end
            for h = obj.hGUI, close(h); end
            obj.pStatus = obj.CONST_STATUS.undefined;
        end
        
        function reload(obj)
            if obj.pStatus < obj.CONST_STATUS.loaded
                addpath(sprintf(['%s' pathsep],obj.toolInPath{:}))
                obj.pStatus = obj.CONST_STATUS.loaded;
            end
            if obj.showGUI, for h = obj.hGUI, set(h,'visible','on'); end; end
        end
        
        function unload(obj)
            if obj.pStatus > obj.CONST_STATUS.unloaded
                warning('%s''s folders (and subfolders) will be removed from the MATLAB path',class(obj));
                rmpath(sprintf(['%s' pathsep],obj.toolInPath{:}))
                obj.pStatus = obj.CONST_STATUS.unloaded;
            end
            if obj.showGUI, for h = obj.hGUI, set(h,'visible','off'); end; end
        end
    end
end