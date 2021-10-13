# Toolboxes

The toolbox in _aa_ is a framework for interacting with extrenal tools. They ensure tight integration and transparent control over their operation. The main concept is that an external tool SHOULD be present in only where and when it is needed to avoid ambiguity due to unintended "shadowing" caused by similar or same function names. This is an issue primarily for MATLAB-based tools; therefore, the toolbox framework also primarily concerns these tools.

## Components of the framework

To integrate an external tool for a particular module, a developer MUST
1. Write an __interface object__
2. Ensure that the corresponing __module__ checks, loads, and unloads the tool.

Also, if a user wants to use this module, they MUST
1. manually download the external tool 
2. Add a corresponding toolbox entry to the site-specific) __parameterset__

### Interface object

The interface is a MATLAB object derived from the [toolboxClass](./toolboxClass.m)

The minimum interface MUST 
- have name of `<(abbreviated) tool name in all lowercase>Class`
- be a subclass of `toolboxClass`
- contain a protected property `hGUI` to store handles to GUI even if the toolbox has none
- contain an initialisation method parsing the arguments and setting argument defaults
- contain a `load` method adding the requires folders to the MATLAB path. This `load` methods overrides and calls `toolboxClass.load`.

For example ([fwsClass](./fwsClass.m)):

```matlab
% aa toolbox interface for Fusion-Watershed (Computational, Cognitive and Clinical Neuroimaging Laboratory, ICL, London, UK)
classdef fwsClass < toolboxClass % fwsClass is a subclass of toolboxClass 
    properties (Access = protected)
        hGUI = [] % protected compuslory store for GUI handles
    end
    
    methods
        function obj = fwsClass(path,varargin)
            defaultAddToPath = false; % default for doAddToPath, typically true only for SPM because its function are used throughout aa
            
            argParse = inputParser;
            argParse.addRequired('path',@ischar);
            argParse.addParameter('name','',@ischar);
            argParse.addParameter('doAddToPath',defaultAddToPath,@(x) islogical(x) || isnumeric(x));
            argParse.parse(path,varargin{:});
            
            % calling overriden method with preset arguments
            %     name                      - Name of the tool corresponding to the name of the object (i.e. 'fws' in this case)
            %     path                      - Path to the tool's main folder
            %     doAddToPath               - Add tool to the path upon initialisation (not the same as loading)
            %     workspaceVariableNames    - List of variables to be stored, typically an empty cell array for no variables
            obj = obj@toolboxClass(argParse.Results.name,argParse.Results.path,argParse.Results.doAddToPath,{}); 
        end
        
        function load(obj)
            addpath(obj.toolPath); % adding tool folder to the MATLAB path
            
            load@toolboxClass(obj) % calling overriden method
        end
    end
end
```

### Module

The corresponding module manages the tool via the interface. The minimum management includes (example from `aamod_tdt_decode`)
- checking the presence of the toolbox during workflow initialisation (`aa_init`)
```matlab
    case 'checkrequirements'
        if ~aas_cache_get(aap,'tdt'), aas_log(aap,true,'TDT is not found'); end
```

- loading the the toolbox at the beginning of execution (`'doit'`)
```matlab
        TDT.load;
```

- unloading the toolbox at the end of the execution (`'doit'`)
```matlab
        TDT.unload;
```

### Parameterset

The site-specific parameterset contains a corresponding toolbox entry for the tool. The mimumum entry includes name (corresponding to the interface name) and path. Multiple entries can be added right after each other. The order does not matter.

For example ([aap_parameters_defaults_UoS](../../aa_parametersets/aap_parameters_defaults_UoS.xml)):

```xml
            <toolbox desc='Toolbox with implemented interface in extrafunctions/toolboxes' ui='custom'>
                <name desc='Name corresponding to the name of the interface without the "Class" suffix' ui='text'>tdt</name>
                <dir ui='dir'>/users/psychology01/software/decoding_toolbox</dir>
            </toolbox>
            <toolbox desc='Toolbox with implemented interface in extrafunctions/toolboxes' ui='custom'>
                <name desc='Name corresponding to the name of the interface without the "Class" suffix' ui='text'>fws</name>
                <dir ui='dir'>/users/psychology01/software/matlab-fws</dir>
            </toolbox>
            <toolbox desc='Toolbox with implemented interface in extrafunctions/toolboxes' ui='custom'>
                <name desc='Name corresponding to the name of the interface without the "Class" suffix' ui='text'>bnv</name>
                <dir ui='dir'>/users/psychology01/software/BrainNetViewer</dir>
            </toolbox>
            <toolbox desc='Toolbox with implemented interface in extrafunctions/toolboxes' ui='custom'>
                <name desc='Name corresponding to the name of the interface without the "Class" suffix' ui='text'>conn</name>
                <dir ui='dir'>/users/psychology01/software/CONN</dir>
            </toolbox>
```

## Optional extras