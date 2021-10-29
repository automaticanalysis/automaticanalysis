classdef connClass < toolboxClass
    properties (Access = protected)
        hGUI = []% GUI handles
        CONN_x0
    end
       
    methods
        function obj = connClass(path,varargin)
            defaultAddToPath = false;
            
            argParse = inputParser;
            argParse.addRequired('path',@ischar);
            argParse.addParameter('name','',@ischar);
            argParse.addParameter('doAddToPath',defaultAddToPath,@(x) islogical(x) || isnumeric(x));
            argParse.parse(path,varargin{:});
            
            vars = {...
                '{"name": "CONN_x", "attributes": ["global"]}'...
                '{"name": "CONN_h", "attributes": ["global"]}'...
                '{"name": "CONN_gui", "attributes": ["global"]}'...
                };
            
            obj = obj@toolboxClass(argParse.Results.name,argParse.Results.path,argParse.Results.doAddToPath,vars);
        end
        
        function load(obj)
            addpath(obj.toolPath);
            if obj.showGUI
                conn;
                obj.hGUI = [CONN_h.screen.hfig CONN_h.screen.hlog];
            else
                conn('init');
            end
            
            global CONN_x
            if ~obj.showGUI, CONN_x.gui = 0; end
            
            % clear
            CONN_x.Setup.rois.names = {};
            CONN_x.Setup.rois.files = {};
            CONN_x.Setup.rois.dimensions = {};
            CONN_x.Setup.rois.mask = [];
            CONN_x.Setup.rois.subjectspecific = [];
            CONN_x.Setup.rois.sessionspecific = [];
            CONN_x.Setup.rois.multiplelabels = [];
            CONN_x.Setup.rois.regresscovariates = [];
            CONN_x.Setup.rois.unsmoothedvolumes = [];
            CONN_x.Setup.rois.weighted = [];
            CONN_x.Setup.conditions.names = {''};
            CONN_x.Analyses(1) = [];
            CONN_x.Results(1) = [];
            
            obj.CONN_x0 = CONN_x;
            
            load@toolboxClass(obj)
        end
        
        function reset(obj)            
            global CONN_x
            CONN_x = obj.CONN_x0;
        end
    end
end