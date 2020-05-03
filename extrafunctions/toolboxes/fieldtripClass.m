classdef fieldtripClass < toolboxClass
    methods
        function obj = fieldtripClass(path,do_add_to_path)
            if nargin < 2, do_add_to_path = false; end
            obj = obj@toolboxClass(path,do_add_to_path);
        end
        
        function init(obj)
            addpath(obj.tool_path);
            ft_defaults
            spmver = '';
            try spmver = spm('ver'); catch, warning('SPM is not detected'); end
            if ~isempty(spmver)
                global ft_default
                ft_default.trackcallinfo = 'no';
                ft_default.showcallinfo = 'no';
                ft_default.spmversion = lower(spmver);
            end
            
            init@toolboxClass(obj)
        end
        
    end
end