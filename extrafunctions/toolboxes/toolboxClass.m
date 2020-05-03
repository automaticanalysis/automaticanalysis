classdef toolboxClass < handle
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
        p_status = -1
        
        tool_path = ''
        tool_in_path = {}
    end
    
    methods
        function obj = toolboxClass(path,do_add_to_path)
            obj.tool_path = path;
            if do_add_to_path, addpath(obj.tool_path); end
            obj.p_status = obj.CONST_STATUS.defined;
        end
        
        function delete(obj)
            obj.close;
        end
        
        function val = get.status(obj)
            for f = fieldnames(obj.CONST_STATUS)'
                if obj.CONST_STATUS.(f{1}) == obj.p_status
                    val = f{1};
                    break;
                end
            end
        end
        
        function init(obj)
            if obj.p_status < obj.CONST_STATUS.loaded
                p = split(path,pathsep);
                obj.tool_in_path = p(cellfun(@(x) contains(x,obj.tool_path), p));
                obj.p_status = obj.CONST_STATUS.loaded;
            end
        end
        
        function close(obj)
            obj.remove_from_path;    
            obj.p_status = obj.CONST_STATUS.undefined;
        end
        
        function add_to_path(obj)
            if obj.p_status < obj.CONST_STATUS.loaded
                addpath(sprintf(['%s' pathsep],obj.tool_in_path{:}))
                obj.p_status = obj.CONST_STATUS.loaded;
            end
        end
        
        function remove_from_path(obj)
            if obj.p_status > obj.CONST_STATUS.unloaded
                rmpath(sprintf(['%s' pathsep],obj.tool_in_path{:}))
                obj.p_status = obj.CONST_STATUS.unloaded;
            end
        end
    end
end