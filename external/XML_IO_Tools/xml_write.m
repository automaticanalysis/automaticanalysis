function DOMnode = xml_write(filename, tree, RootName, Pref)
%XML_WRITE  Writes Matlab data structures to XML file
%
% DESCRIPTION
% xml_write( filename, tree) Converts Matlab data structure 'tree' containing
% cells, structs, numbers and strings to Document Object Model (DOM) node
% tree, then saves it to XML file 'filename' using Matlab's xmlwrite
% function. Optionally one can also use alternative version of xmlwrite
% function which directly calls JAVA functions for XML writing without
% MATLAB middleware. This function is provided as a patch to existing
% bugs in xmlwrite (in R2006b).
%
% xml_write(filename, tree, RootName, Pref) allows you to specify
% additional preferences about file format
%
% DOMnode = xml_write([], tree) same as above except that DOM node is
% not saved to the file but returned.
%
% INPUT
%   filename     file name
%   tree         Matlab structure tree to store in xml file.
%   RootName     String with XML tag name used for root (top level) node
%                Optionally it can be a string cell array storing: Name of
%                root node, document "Processing Instructions" data and
%                document "comment" string
%   Pref         Other preferences:
%     Pref.ItemName - default 'item' -  name of a special tag used to
%                     itemize cell arrays
%     Pref.XmlEngine - let you choose the XML engine. Currently default is
%       'Xerces', which is using directly the apache xerces java file.
%       Other option is 'Matlab' which uses MATLAB's xmlwrite and its
%       XMLUtils java file. Both options create identical results except in
%       case of CDATA sections where xmlwrite fails.
%     Pref.CellItem - default 'true' - allow cell arrays to use 'item'
%       notation. See below.
%    Pref.RootOnly - default true - output variable 'tree' corresponds to
%       xml file root element, otherwise it correspond to the whole file.
%     Pref.StructItem - default 'true' - allow arrays of structs to use
%       'item' notation. For example "Pref.StructItem = true" gives:
%         <a>
%           <b>
%             <item> ... <\item>
%             <item> ... <\item>
%           <\b>
%         <\a>
%       while "Pref.StructItem = false" gives:
%         <a>
%           <b> ... <\b>
%           <b> ... <\b>
%         <\a>
%
%
% Several special xml node types can be created if special tags are used
% for field names of 'tree' nodes:
%  - node.CONTENT - stores data section of the node if other fields
%    (usually ATTRIBUTE are present. Usually data section is stored
%    directly in 'node'.
%  - node.ATTRIBUTE.name - stores node's attribute called 'name'.
%  - node.COMMENT - create comment child node from the string. For global
%    comments see "RootName" input variable.
%  - node.PROCESSING_INSTRUCTIONS - create "processing instruction" child
%    node from the string. For global "processing instructions" see
%    "RootName" input variable.
%  - node.CDATA_SECTION - stores node's CDATA section (string). Only works
%    if Pref.XmlEngine='Xerces'. For more info, see comments of F_xmlwrite.
%  - other special node types like: document fragment nodes, document type
%    nodes, entity nodes and notation nodes are not being handled by
%    'xml_write' at the moment.
%
% OUTPUT
%   DOMnode      Document Object Model (DOM) node tree in the format
%                required as input to xmlwrite. (optional)
%
% EXAMPLES:
%   MyTree=[];
%   MyTree.MyNumber = 13;
%   MyTree.MyString = 'Hello World';
%   xml_write('test.xml', MyTree);
%   type('test.xml')
%   %See also xml_tutorial.m
%
% See also
%   xml_read, xmlread, xmlwrite
%
% Written by Jarek Tuszynski, SAIC, jaroslaw.w.tuszynski_at_saic.com

%% Check Matlab Version
v = ver('MATLAB');
v = str2double(regexp(v.Version, '\d.\d','match','once'));
if (v<7)
  error('Your MATLAB version is too old. You need version 7.0 or newer.');
end

%% default preferences
DPref.ItemName  = 'item'; % name of a special tag used to itemize cell arrays
DPref.StructItem = true;  % allow arrays of structs to use 'item' notation
DPref.CellItem   = true;  % allow cell arrays to use 'item' notation
DPref.XmlEngine  = 'Matlab';  % use matlab provided XMLUtils
%DPref.XmlEngine  = 'Xerces';  % use Xerces xml generator directly
RootOnly         = true;  % Input is root node only
GlobalProcInst = [];
GlobalComment  = [];
GlobalDocType  = [];

%% read user preferences
if (nargin>3)
  if (isfield(Pref, 'ItemName'  )), DPref.ItemName   = Pref.ItemName;   end
  if (isfield(Pref, 'StructItem')), DPref.StructItem = Pref.StructItem; end
  if (isfield(Pref, 'CellItem'  )), DPref.CellItem   = Pref.CellItem;   end
  if (isfield(Pref, 'XmlEngine' )), DPref.XmlEngine  = Pref.XmlEngine;  end
  if (isfield(Pref, 'RootOnly'  )), RootOnly         = Pref.RootOnly;   end
end
if (nargin<3 || isempty(RootName)), RootName=inputname(2); end
if (isempty(RootName)), RootName='ROOT'; end
if (iscell(RootName)) % RootName also stores global text node data
  rName = RootName;
  RootName = char(rName{1});
  if (length(rName)>1), GlobalProcInst = char(rName{2}); end
  if (length(rName)>2), GlobalComment  = char(rName{3}); end
  if (length(rName)>3), GlobalDocType  = char(rName{4}); end
end
if(~RootOnly && isstruct(tree))  % if struct than deal with each field separatly
  fields = fieldnames(tree);
  for i=1:length(fields)
    field = fields{i};
    x = tree(1).(field);
    if (strcmp(field, 'COMMENT'))
      GlobalComment = x;
    elseif (strcmp(field, 'PROCESSING_INSTRUCTION'))
      GlobalProcInst = x;
    elseif (strcmp(field, 'DOCUMENT_TYPE'))
      GlobalDocType = x;
    else
      RootName = field;
      t = x;
    end
  end
  tree = t;
end

%% Initialize jave object that will store xml data structure
RootName = varName2str(RootName);
if (~isempty(GlobalDocType))
  %   n = strfind(GlobalDocType, ' ');
  %   if (~isempty(n))
  %     dtype = com.mathworks.xml.XMLUtils.createDocumentType(GlobalDocType);
  %   end
  %   DOMnode = com.mathworks.xml.XMLUtils.createDocument(RootName, dtype);
  warning('xml_io_tools:write:docType', ...
   'DOCUMENT_TYPE node was encountered which is not supported yet. Ignoring.');
end
DOMnode = com.mathworks.xml.XMLUtils.createDocument(RootName);


%% Use recursive function to convert matlab data structure to XML
root = DOMnode.getDocumentElement;
struct2DOMnode(DOMnode, root, tree, DPref.ItemName, DPref);

%% Remove the only child of the root node
root   = DOMnode.getDocumentElement;
Child  = root.getChildNodes; % create array of children nodes
nChild = Child.getLength;    % number of children
if (nChild==1)
  node = root.removeChild(root.getFirstChild);
  while(node.hasChildNodes)
    root.appendChild(node.removeChild(node.getFirstChild));
  end
  while(node.hasAttributes)            % copy all attributes
    root.setAttributeNode(node.removeAttributeNode(node.getAttributes.item(0)));
  end
end

%% Save exotic Global nodes
if (~isempty(GlobalComment))
  DOMnode.insertBefore(DOMnode.createComment(GlobalComment), DOMnode.getFirstChild());
end
if (~isempty(GlobalProcInst))
  n = strfind(GlobalProcInst, ' ');
  if (~isempty(n))
    proc = DOMnode.createProcessingInstruction(GlobalProcInst(1:(n(1)-1)),...
      GlobalProcInst((n(1)+1):end));
    DOMnode.insertBefore(proc, DOMnode.getFirstChild());
  end
end
% if (~isempty(GlobalDocType))
%   n = strfind(GlobalDocType, ' ');
%   if (~isempty(n))
%     dtype = DOMnode.createDocumentType(GlobalDocType);
%     DOMnode.insertBefore(dtype, DOMnode.getFirstChild());
%   end
% end

%% save java DOM tree to XML file
if (~isempty(filename))
  if (strcmpi(DPref.XmlEngine, 'Xerces'))
    xmlwrite_xerces(filename, DOMnode);
  else
    xmlwrite(filename, DOMnode);
  end
end


%% =======================================================================
%  === struct2DOMnode Function ===========================================
%  =======================================================================
function [] = struct2DOMnode(xml, parent, s, name, Pref)
% struct2DOMnode is a recursive function that converts matlab's structs to
% DOM nodes.
% INPUTS:
%  xml - jave object that will store xml data structure
%  parent - parent DOM Element
%  s - Matlab data structure to save
%  name - name to be used in xml tags describing 's'
%  Pref - preferenced
name = varName2str(name);
ItemName = Pref.ItemName;
if (ischar(s) && min(size(s))>1) % if 2D array of characters
  s=cellstr(s);                  % than convert to cell array
end
while (iscell(s) && length(s)==1), s = s{1}; end
nItem = length(s);
if (iscell(s)) % if this is a cell array
  if (nItem==1 || strcmp(name, 'CONTENT') || ~Pref.CellItem)
    if (strcmp(name, 'CONTENT')), CellName = ItemName;  % use 'item' notation <item> ... <\item>
    else                          CellName = name; end  % don't use 'item' notation <a> ... <\a>
    for iItem=1:nItem   % save each cell separatly
      struct2DOMnode(xml, parent, s{iItem}, CellName, Pref); % recursive call
    end
  else % use 'item' notation  <a> <item> ... <\item> <\a>
    node = xml.createElement(name);
    for iItem=1:nItem   % save each cell separatly
      struct2DOMnode(xml, node, s{iItem}, ItemName , Pref); % recursive call
    end
    parent.appendChild(node);
  end
elseif (isstruct(s))  % if struct than deal with each field separatly
  fields = fieldnames(s);
  % if array of structs with no attributes than use 'items' notation
  if (nItem>1 && Pref.StructItem && ~isfield(s,'ATTRIBUTE') )
    node = xml.createElement(name);
    for iItem=1:nItem
      struct2DOMnode(xml, node, s(iItem), ItemName, Pref ); % recursive call
    end
    parent.appendChild(node);
  else % otherwise save each struct separatelly
    for j=1:nItem
      node = xml.createElement(name);
      for i=1:length(fields)
        field = fields{i};
        x = s(j).(field);
        %if (isempty(x)), continue; end
        if (iscell(x) && (strcmp(field, 'COMMENT') || ...
            strcmp(field, 'CDATA_SECTION') || ...
            strcmp(field, 'PROCESSING_INSTRUCTION')))
          for k=1:length(x) % if nodes that should have strings have cellstrings
            struct2DOMnode(xml, node, x{k}, field, Pref ); % recursive call will modify 'node'
          end
        elseif (strcmp(field, 'ATTRIBUTE')) % set attributes of the node
          if (~isstruct(x))
            warning('xml_io_tools:write:badAttribute', ...
              'Struct field named ATTRIBUTE encountered which was not a struct. Ignoring.');
            continue;
          end
          attName = fieldnames(x);       % get names of all the attributes
          for k=1:length(attName)        % attach them to the node
            att = xml.createAttribute(varName2str(attName(k)));
            att.setValue(var2str(x.(attName{k})));
            node.setAttributeNode(att);
          end
        else                             % set children of the node
          struct2DOMnode(xml, node, x, field, Pref ); % recursive call will modify 'node'
        end
      end  % end for i=1:nFields
      parent.appendChild(node);
    end  % end for j=1:nItem
  end
else  % if not a struct and not a cell than it is a leaf node
  if (strcmp(name, 'CONTENT'))
    txt = xml.createTextNode(var2str(s)); % ... than it can be converted to text
    parent.appendChild(txt);
  elseif (strcmp(name, 'COMMENT'))   % create comment node
    if (ischar(s))
      com = xml.createComment(s);
      parent.appendChild(com);
    else
      warning('xml_io_tools:write:badComment', ...
        'Struct field named COMMENT encountered which was not a string. Ignoring.');
    end
  elseif (strcmp(name, 'CDATA_SECTION'))   % create CDATA Section
    if (ischar(s))
      cdt = xml.createCDATASection(s);
      parent.appendChild(cdt);
    else
      warning('xml_io_tools:write:badCData', ...
        'Struct field named CDATA_SECTION encountered which was not a string. Ignoring.');
    end
  elseif (strcmp(name, 'PROCESSING_INSTRUCTION')) % set attributes of the node
    OK = false;
    if (ischar(s))
      n = strfind(s, ' ');
      if (~isempty(n))
        proc = xml.createProcessingInstruction(s(1:(n(1)-1)),s((n(1)+1):end));
        parent.insertBefore(proc, parent.getFirstChild());
        OK = true;
      end
    end
    if (~OK)
      warning('xml_io_tools:write:badProcInst', ...
        ['Struct field named PROCESSING_INSTRUCTION need to be',...
        ' a string, for example: xml-stylesheet type="text/css" ', ...
        'href="myStyleSheet.css". Ignoring.']);
    end
  else % I guess it is a regular text leaf node
    txt  = xml.createTextNode(var2str(s));
    node = xml.createElement(name);
    node.appendChild(txt);
    parent.appendChild(node);
  end
end

%% =======================================================================
%  === var2str Function ==================================================
%  =======================================================================
function str = var2str(s)
% convert matlab variables to a sting
if (isnumeric(s) || islogical(s))
  dim = size(s);
  if (min(dim)<=1 || length(dim)>2) % if 1D or 3D array
    s=s(:); s=s.';            % convert to 1D array
    str=num2str(s);           % convert array of numbers to string
  else                        % if a 2D array
    s=mat2str(s);             % convert matrix to a string
    str=regexprep(s,';',';\n');
  end
elseif iscell(s)
  str = char(s{1});
  for i=2:length(s)
    str = [str, ' ', char(s{i})];
  end
else
  str = char(s);
end
str=str(:); str=str.';            % make sure this is a row vector of char's
if (length(str)>1)
  str(str<32|str==127)=' ';       % convert no-printable characters to spaces
  str = strtrim(str);             % remove spaces from begining and the end
  str = regexprep(str,'\s+',' '); % remove multiple spaces
end

%% =======================================================================
%  === var2Namestr Function ==============================================
%  =======================================================================
function str = varName2str(str)
% convert matlab variable names to a sting
str = char(str);
p   = strfind(str,'0x');
if (~isempty(p))
  for i=1:length(p)
    before = str( p(i)+(0:3) );          % string to replace
    after  = char(hex2dec(before(3:4))); % string to replace with
    str = regexprep(str,before,after, 'once', 'ignorecase');
    p=p-3; % since 4 characters were replaced with one - compensate
  end
end
str = regexprep(str,'_COLON_',':', 'once', 'ignorecase');
str = regexprep(str,'_DASH_' ,'-', 'once', 'ignorecase');

