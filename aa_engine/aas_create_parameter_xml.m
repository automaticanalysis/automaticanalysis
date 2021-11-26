function parameter_xml_filename = aas_create_parameter_xml(parameter_xml_filename, useGUI, varargin)
% Seed a new parameters file
% parameter_xml_filename = aas_create_parameter_xml(parameter_xml_filename, useGUI, varargin)
%
% Through user interaction the file and directory names are acquired to be
% used for and in the parameter xml file.
%
% Inputs:
%  - parameter_xml_filename: Base file name for the new parameter xml file.
% Can be changed during this function when the user is asked where to store
% the file. Input value will be ignored when use_default_filename = true
%  - useGUI: Whether to use graphical dialogs to interact with the user.
% Command window interaction will be used when set to false.
%  - varargin: Parameter-value pairs
%  -- use_default_location (default false): When true, store the new
%  parameter file in the default (as specified by aa) user config folder
%  -- use_default_filename (default false): When true, use the default user
%  parameter file name (as specified by aa) for the new parameter file

%% Parse inputs
p = inputParser;
addParameter(p, 'use_default_location', false, @islogical);
addParameter(p, 'use_default_filename', false, @islogical);
parse(p,varargin{:});
use_default_location = p.Results.use_default_location;
use_default_filename = p.Results.use_default_filename;

aa_info = aaClass('nopath', 'nogreet');
if use_default_filename
    parameter_xml_filename = aa_info.parameter_xml_filename;
end

%% Which parameter set to use as seed
[aahome,~,~] = fileparts(which('aa_ver5'));
defaultdir = fullfile(aahome,'aa_parametersets');
ui_msg = 'Select parameter set that will be used as seed';
[seedparam, rootpath] = userinput(...
    'uigetfile',{'*.xml','All Parameter Files' },ui_msg,defaultdir,'GUI',useGUI);
assert(ischar(seedparam), 'Exiting, user cancelled');
is_base_defaults_file = strcmp(seedparam, 'aap_parameters_defaults.xml');
seedparam = fullfile(rootpath, seedparam);

%% Where to store the new parameters file
destination = fullfile(aa_info.configdir, parameter_xml_filename);
if ~use_default_location
    ui_msg = 'Location where the parameters file will be saved';
    [parameter_xml_filename, rootpath] = userinput(...
        'uiputfile',{'*.xml','All Parameter Files' },ui_msg,destination,'GUI',useGUI);
    assert(ischar(parameter_xml_filename), 'Exiting, user cancelled');
    destination = fullfile(rootpath, parameter_xml_filename);
end

%% Get value for acq_details.root
% Initialise the save dialogue in the current aap.acq_details.root if specified
xml = xml_read(seedparam);
analysisroot = aas_expandpathbyvars(xml.acq_details.root.CONTENT);
previous = '';
while ~isempty(analysisroot) && ~strcmp(previous, analysisroot)
    if exist(analysisroot, "dir")
        break
    end
    previous = analysisroot;
    analysisroot = fileparts(analysisroot);
end
ui_msg = 'Location where intermediate and final analysis results will be stored';
analysisroot = userinput('uigetdir',analysisroot,ui_msg,'GUI',useGUI);
assert(ischar(analysisroot), 'Exiting, user cancelled');

%% Get values for other required directories.
rawdataroot = '';
spmroot = '';
if is_base_defaults_file
    % Only when the base file was selected as seed. When a specific
    % location's file was selected, assume that these values are already ok
    % for that location.
    ui_msg = 'Directory where raw input data can be found / will be stored';
    rawdataroot = userinput('uigetdir',analysisroot,ui_msg,'GUI',useGUI);
    assert(ischar(rawdataroot), 'Exiting, user cancelled');

    ui_msg = 'Root directory of spm installation';
    spmroot = userinput('uigetdir',analysisroot,ui_msg,'GUI',useGUI);
    assert(ischar(spmroot), 'Exiting, user cancelled');
end

%% Generate new parameters file
create_minimalXML(seedparam, destination, analysisroot, rawdataroot, spmroot);
assert(exist(destination,'file')>0,'Failed to create %s', destination);

%% Final check and messaging
% The file should now be on the path. But check, it might not be e.g. if
% aa was not added to the path properly before calling this function.
if ~exist(parameter_xml_filename,'file')
    msg = sprintf('Could not find %s - Are you sure it is on your path?', parameter_xml_filename);
    ME = MException('aas_create_parameter_xml:notOnPath', msg);
    throw(ME)
end

msg = sprintf('New parameter set in %s has been created.\nYou may need to edit this file further to reflect local configuration.',destination);
if useGUI
    h = msgbox(msg,'New parameters file','Warn');
    waitfor(h);
else
    fprintf('\n%s\n',msg);
end

end

%% create_minimalXML
function create_minimalXML(seedparam,destination,analysisroot, rawdataroot, spmroot)

docNode = com.mathworks.xml.XMLUtils.createDocument('aap');
aap = docNode.getDocumentElement;
aap.setAttribute('xmlns:xi','http://www.w3.org/2001/XInclude');

seed = docNode.createElement('xi:include');
seed.setAttribute('href',seedparam);
seed.setAttribute('parse','xml');
aap.appendChild(seed);

local = docNode.createElement('local');
aap.appendChild(local);

if ~isempty(rawdataroot) ||  ~isempty(spmroot)
    directory_conventions = docNode.createElement('directory_conventions');
    local.appendChild(directory_conventions);

    if ~isempty(rawdataroot)
        rawdatadir = docNode.createElement('rawdatadir');
        rawdatadir.setAttribute('desc','Root on local machine for processed data');
        rawdatadir.setAttribute('ui','dir');
        rawdatadir.appendChild(docNode.createTextNode(rawdataroot));
        directory_conventions.appendChild(rawdatadir);
    end

    if ~isempty(spmroot)
        toolbox = docNode.createElement('toolbox');
        toolbox.setAttribute('desc','Toolbox with implemented interface in extrafunctions/toolboxes');
        toolbox.setAttribute('ui','custom');
        directory_conventions.appendChild(toolbox);

        spmname = docNode.createElement('name');
        spmname.setAttribute('desc','Name corresponding to the name of the interface without the "Class" suffix');
        spmname.setAttribute('ui','text');
        spmname.appendChild(docNode.createTextNode('spm'));
        toolbox.appendChild(spmname);

        spmdir = docNode.createElement('dir');
        spmdir.setAttribute('ui','dir_list');
        spmdir.appendChild(docNode.createTextNode(spmroot));
        toolbox.appendChild(spmdir);
    end
end

acq_details = docNode.createElement('acq_details');
local.appendChild(acq_details);

root = docNode.createElement('root');
root.setAttribute('desc','Root on local machine for processed data');
root.setAttribute('ui','dir');
root.appendChild(docNode.createTextNode(analysisroot));
acq_details.appendChild(root);

xmlwrite(destination,docNode);

end
