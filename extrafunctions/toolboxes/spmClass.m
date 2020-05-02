classdef spmClass < toolboxClass
    methods
        function obj = spmClass(path,do_add_to_path)
            if nargin < 2, do_add_to_path = false; end
            obj = obj@toolboxClass(path,do_add_to_path);
        end
        
        function init(obj)
            addpath(obj.tool_path);
            spm_jobman('initcfg');
            
            init@toolboxClass(obj)
        end
        
    end
end