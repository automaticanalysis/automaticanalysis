% N.B.: This is not an integration of FreeSurfer, yet.
% This class, as is, integrates labels (name and MNI coord) for sereval
% FreeSurfer atlases. These labels are bundled with aa as fs_mods.
classdef freesurferClass < toolboxClass
    properties
        labels
    end
    
    properties (Access = protected)
        hGUI = []% GUI handles
    end
    
    methods
        function obj = freesurferClass(path,varargin)
            defaultAddToPath = false;
            
            argParse = inputParser;
            argParse.addRequired('path',@ischar);
            argParse.addParameter('name','',@ischar);
            argParse.addParameter('doAddToPath',defaultAddToPath,@(x) islogical(x) || isnumeric(x));
            argParse.parse(path,varargin{:});
            
            obj = obj@toolboxClass(argParse.Results.name,argParse.Results.path,argParse.Results.doAddToPath,{});
        end
        
        function load(obj)
            load@toolboxClass(obj)
            
            labelDir = fullfile(fileparts(mfilename('fullpath')),'freesurfer_mods','labels');
            for fn = cellstr(spm_select('FPList',labelDir,'^aparc_.*'))'
                atlasName = strrep(spm_file(fn{1},'basename'),'aparc_','');
                T = readtable(fn{1},'FileType', 'text');
                T.Properties.VariableNames = {'Structure' 'X' 'Y' 'Z'};
                T.Properties.VariableUnits = {'' 'mm' 'mm' 'mm'};
                obj.labels.(atlasName) = T;
            end
        end
    end
end