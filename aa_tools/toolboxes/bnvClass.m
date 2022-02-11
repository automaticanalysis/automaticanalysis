classdef bnvClass < toolboxClass
    properties (Access = protected)
        hGUI = []% GUI handles
    end
    
    properties
        displaySettings = struct(...
            'view','LMD',...
            'edgeSize', 0.5,...
            'LUT','hot'...
            )
    end
    
    properties (Dependent)
        BNVSettings
    end
    
    methods
        function obj = bnvClass(path,varargin)
            defaultAddToPath = false;
            
            argParse = inputParser;
            argParse.addRequired('path',@ischar);
            argParse.addParameter('name','',@ischar);
            argParse.addParameter('doAddToPath',defaultAddToPath,@(x) islogical(x) || isnumeric(x));
            argParse.parse(path,varargin{:});
            
            obj = obj@toolboxClass(argParse.Results.name,argParse.Results.path,argParse.Results.doAddToPath,{});
        end
        
        function load(obj)
            addpath(obj.toolPath);
            addpath(fullfile(obj.toolPath,'BrainNet_spm12_files'));

            load@toolboxClass(obj)
        end
        
        function out = get.BNVSettings(obj)
            VIEW = {...
                'S' % single
                'F' % full
                'LM' % lateral and medial view
                'LMV' % lateral, medial and ventral view
                'LMD' % lateral, medial and dorsal view
                };
            
            out = struct;
            out.lot = struct;
            out.lot.view = find(strcmp(VIEW,obj.displaySettings.view));
            out.edg = struct;
            out.edg.size = 1;
            out.edg.size_size = obj.displaySettings.edgeSize;
            out.edg.color = 2;
            out.edg.color_map = 3;
            cm = feval(obj.displaySettings.LUT,128); % avoid edges with too light or dark colours
            out.edg.CM = cm(33:96,:);
            out.edg.directed = 1;
        end
    end
end