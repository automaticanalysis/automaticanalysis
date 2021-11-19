function parameter_xml_filename = aas_create_parameter_xml(parameter_xml_filename, useGUI)
% if we made it here, we are seeding a new parameters file
% we have default parameters
defaultdir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'aa_parametersets');
[seedparam, rootpath] = userinput('uigetfile',{'*.xml','All Parameter Files' },'Desired seed parameter',defaultdir,'GUI',useGUI);
assert(ischar(seedparam), 'exiting');
seedparam = fullfile(rootpath, seedparam);

% initialise the save dialogue in the current aap.acq_details.root if specified
xml=xml_read(seedparam);
configdir = fullfile(getenv('HOME'),'.aa');
% generate new parameters file N.B.: in networks with shared resources
% average user may not be able to write into aa_paremetersets
[parameter_xml_filename, rootpath] = userinput('uiputfile',{'*.xml','All Parameter Files' },...
    'Location of the parameters file',fullfile(configdir, parameter_xml_filename),'GUI',useGUI);
assert(ischar(parameter_xml_filename), 'exiting');
destination = fullfile(rootpath, parameter_xml_filename);

analysisroot = aas_expandpathbyvars(xml.acq_details.root.CONTENT);
% TODO: Why the makedir here, if we are going to ask the user for a
% location? Now creates a (chain-of-)directory, that then might not be
% used.
% Only use this if available, else use the first directory that is. If
% nothing, use fullfile(getenv('HOME'),'.aa').
aas_makedir([], analysisroot);
analysisroot = userinput('uigetdir',analysisroot,'Location of analyses by default','GUI',useGUI); % TODO: Missing an assert on the returned value

create_minimalXML(seedparam, destination, analysisroot);
assert(exist(destination,'file')>0,'failed to create %s', parameter_xml_filename);

% N.B. we don't actually modify parameter_xml_filename - it should now be on the path. But
% let's double check. It might not be e.g. if you haven't actually added AA to your
% path properly before calling this function.
assert(exist(parameter_xml_filename,'file')>0, ...
    'could not find %s - Are you sure it is on your path?',...
    parameter_xml_filename);

if useGUI
    h = msgbox(sprintf('New parameter set in %s has been created.\nYou may need to edit this file further to reflect local configuration.',destination),'New parameters file','Warn');
    waitfor(h);
else
    fprintf('\nNew parameter set in %s has been created.\nYou may need to edit this file further to reflect local configuration.\n',destination);
end

end

%% create_minimalXML
function create_minimalXML(seedparam,destination,analysisroot)

if nargin < 3, analysisroot = fileparts(destination); end

docNode = com.mathworks.xml.XMLUtils.createDocument('aap');
aap = docNode.getDocumentElement;
aap.setAttribute('xmlns:xi','http://www.w3.org/2001/XInclude');

seed = docNode.createElement('xi:include');
seed.setAttribute('href',seedparam);
seed.setAttribute('parse','xml');
aap.appendChild(seed);

local = docNode.createElement('local');
aap.appendChild(local);

acq_details = docNode.createElement('acq_details');
local.appendChild(acq_details);

root = docNode.createElement('root');
root.setAttribute('desc','Root on local machine for processed data');
root.setAttribute('ui','dir');
root.appendChild(docNode.createTextNode(analysisroot));
acq_details.appendChild(root);

xmlwrite(destination,docNode);

end
