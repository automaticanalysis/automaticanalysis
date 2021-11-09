% this requires AROMA-ICA to be installed and added as an aa toolbox
%
% The easiest way to install is prolly:
%
%   % cd /users/abcd1234/tools
%   % sudo git clone https://github.com/maartenmennes/ICA-AROMA.git
%
%   (assuming that the repo link is still valid)
%
% add the correspodning entry to your parameterset
%   <toolbox desc='Toolbox with implemented interface in extrafunctions/toolboxes' ui='custom'>
%       <name desc='Name corresponding to the name of the interface without the "Class" suffix' ui='text'>aroma</name>
%           <dir ui='dir'>/users/abcd1234/tools/ICA-AROMA</dir>
%   </toolbox>
classdef aromaClass < toolboxClass
    properties (Access = protected)
        hGUI = []% GUI handles
    end
    
    methods
        function obj = aromaClass(path,varargin)
            defaultAddToPath = false;
            defaultKeepInPath = false;
            
            argParse = inputParser;
            argParse.addRequired('path',@ischar);
            argParse.addParameter('name','',@ischar);
            argParse.addParameter('doAddToPath',defaultAddToPath,@(x) islogical(x) || isnumeric(x));
            argParse.addParameter('doKeepInPath',defaultKeepInPath,@(x) islogical(x) || isnumeric(x));
            argParse.parse(path,varargin{:});
            
            obj = obj@toolboxClass(argParse.Results.name,argParse.Results.path,argParse.Results.doAddToPath,argParse.Results.doKeepInPath,{});
        end
        
        function load(obj)
            load@toolboxClass(obj)

        end
    end
end