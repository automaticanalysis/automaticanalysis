function [tree, RootName, DOMnode] = xml_read(xmlfile, Pref)
%XML_READ reads xml files and converts them into Matlab's struct tree.
%
% DESCRIPTION
% tree = xml_read(xmlfile) reads 'xmlfile' into data structure 'tree'
%
% tree = xml_read(xmlfile, Pref) reads 'xmlfile' into data structure 'tree'
% according to your preferences
%
% [tree, RootName, DOMnode] = xml_read(xmlfile) get additional information
% about XML file
%
% INPUT:
%  xmlfile	URL or filename of xml file to read
%  Pref     Preferences:
%    Pref.ItemName - default 'item' - name of a special tag used to itemize
%                    cell arrays
%    Pref.ReadAttr - default true - allow reading attributes
%    Pref.ReadSpec - default true - allow reading special nodes
%    Pref.Str2Num  - default true - convert strings that look like numbers
%                   to numbers
%    Pref.NoCells  - default true - force output to have no cell arrays
%    Pref.Debug    - default false - show mode specific error messages
%    Pref.NumLevels- default infinity - how many recursive levels are
%      allowed. Can be used to speed up the function by prunning the tree.
%    Pref.RootOnly - default true - output variable 'tree' corresponds to
%      xml file root element, otherwise it correspond to the whole file.
% OUTPUT:
%  tree         tree of structs and/or cell arrays corresponding to xml file
%  RootName     XML tag name used for root (top level) node.
%               Optionally it can be a string cell array storing: Name of
%               root node, document "Processing Instructions" data and
%               document "comment" string
%  DOMnode      output of xmlread
%
% DETAILS:
% Function xml_read first calls MATLAB's xmlread function and than
% converts its output ('Document Object Model' tree of Java objects)
% to tree of MATLAB struct's. The output is in format of nested structs
% and cells. In the output data structure field names are based on
% XML tags, except in cases when tags produce illegal variable names.
%
% Several special xml node types result in special tags for fields of
% 'tree' nodes:
%  - node.CONTENT - stores data section of the node if other fields are
%    present. Usually data section is stored directly in 'node'.
%  - node.ATTRIBUTE.name - stores node's attribute called 'name'.
%  - node.COMMENT - stores node's comment section (string). For global
%    comments see "RootName" output variable.
%  - node.CDATA_SECTION - stores node's CDATA section (string).
%  - node.PROCESSING_INSTRUCTIONS - stores "processing instruction" child
%    node. For global "processing instructions" see "RootName" output variable.
%  - other special node types like: document fragment nodes, document type
%   nodes, entity nodes, notation nodes and processing instruction nodes
%   will be treated like regular nodes
%
% EXAMPLES:
%   MyTree=[];
%   MyTree.MyNumber = 13;
%   MyTree.MyString = 'Hello World';
%   xml_write('test.xml', MyTree);
%   [tree treeName] = xml_read ('test.xml');
%   disp(treeName)
%   gen_object_display()
%   % See also xml_examples.m
%
% See also:
%   xml_write, xmlread, xmlwrite
%
% Written by Jarek Tuszynski, SAIC, jaroslaw.w.tuszynski_at_saic.com
% References:
%  - Function inspired by Example 3 found in xmlread function.
%  - Output data structures inspired by xml_toolbox structures.
%
% RC 13/2/2010 - took out date conversion option in str2var as this
% erroneously leaves 3-double vectors as strings
%

%% default preferences
DPref.ItemName  = 'item'; % name of a special tag used to itemize cell arrays
DPref.ReadAttr  = true;   % allow reading attributes
DPref.ReadSpec  = true;   % allow reading special nodes: comments, CData, etc.
DPref.Str2Num   = true;   % convert strings that look like numbers to numbers
DPref.NoCells   = true;   % force output to have no cell arrays
DPref.NumLevels = 1e10;   % number of recurence levels
RootOnly        = true;   % return root node  with no top level special nodes
Debug           = false;  % show specific errors (true) or general (false)?
tree            = [];
RootName        = [];

%% Check Matlab Version
v = ver('MATLAB');
v = str2double(regexp(v.Version, '\d.\d','match','once'));
if (v<7.1)
  error('Your MATLAB version is too old. You need version 7.1 or newer.');
end

%% read user preferences
if (nargin>1)
  if (isfield(Pref, 'ItemName' )), DPref.ItemName  = Pref.ItemName;  end
  if (isfield(Pref, 'Str2Num'  )), DPref.Str2Num   = Pref.Str2Num ;  end
  if (isfield(Pref, 'NoCells'  )), DPref.NoCells   = Pref.NoCells ;  end
  if (isfield(Pref, 'NumLevels')), DPref.NumLevels = Pref.NumLevels; end
  if (isfield(Pref, 'ReadAttr' )), DPref.ReadAttr  = Pref.ReadAttr;  end
  if (isfield(Pref, 'ReadSpec' )), DPref.ReadSpec  = Pref.ReadSpec;  end
  if (isfield(Pref, 'RootOnly' )), RootOnly        = Pref.RootOnly;  end
  if (isfield(Pref, 'Debug'    )), Debug           = Pref.Debug   ;  end
end

%% read xml file using Matlab function
parserFactory = javaMethod('newInstance',...
    'javax.xml.parsers.DocumentBuilderFactory');
javaMethod('setXIncludeAware',parserFactory,true);
javaMethod('setNamespaceAware',parserFactory,true);
p = javaMethod('newDocumentBuilder',parserFactory);

if (ischar(xmlfile)) % if xmlfile is a string
  if (Debug)
    DOMnode = xmlread(xmlfile,p);
  else
    try
      DOMnode = xmlread(xmlfile,p);
    catch
      error('Failed to read XML file %s.',xmlfile);
    end
  end
  Node = DOMnode.getFirstChild;
else %if xmlfile is not a string than maybe it is a DOMnode already
  try
    Node = xmlfile.getFirstChild;
    DOMnode = xmlfile;
  catch
    error('Input variable xmlfile has to be a string or DOM node.');
  end
end

%% Find the Root node. Also store data from Global Comment and Processing
%  Instruction nodes, if any.
GlobalTextNodes = cell(1,3);
GlobalProcInst  = [];
GlobalComment   = [];
GlobalDocType   = [];
while (~isempty(Node))
  if (Node.getNodeType==Node.ELEMENT_NODE)
    RootNode=Node;
  elseif (Node.getNodeType==Node.PROCESSING_INSTRUCTION_NODE)
    data   = strtrim(char(Node.getData));
    target = strtrim(char(Node.getTarget));
    GlobalProcInst = [target, ' ', data];
    GlobalTextNodes{2} = GlobalProcInst;
  elseif (Node.getNodeType==Node.COMMENT_NODE)
    GlobalComment = strtrim(char(Node.getData));
    GlobalTextNodes{3} = GlobalComment;
    %   elseif (Node.getNodeType==Node.DOCUMENT_TYPE_NODE)
    %     GlobalTextNodes{4} = GlobalDocType;
  end
  Node = Node.getNextSibling;
end

%% parse xml file through calls to recursive DOMnode2struct function
if (Debug)   % in debuging mode allow crashes
  [tree RootName] = DOMnode2struct(RootNode, DPref, 1);
else         % in normal mode do not allow crashes
  try
    [tree RootName] = DOMnode2struct(RootNode, DPref, 1);
  catch
    error('Unable to parse XML file %s.',xmlfile);
  end
end

%% If there were any Global Text nodes than return them
if (~RootOnly)
  if (~isempty(GlobalProcInst) && DPref.ReadSpec)
    t.PROCESSING_INSTRUCTION = GlobalProcInst;
  end
  if (~isempty(GlobalComment) && DPref.ReadSpec)
    t.COMMENT = GlobalComment;
  end
  if (~isempty(GlobalDocType) && DPref.ReadSpec)
    t.DOCUMENT_TYPE = GlobalDocType;
  end
  t.(RootName) = tree;
  tree=t;
end
if (~isempty(GlobalTextNodes))
  GlobalTextNodes{1} = RootName;
  RootName = GlobalTextNodes;
end

if isfield(tree,'aap') % XInclude used
    tree = mergeStructs(tree.aap,tree.mod);
end

%% =======================================================================
%  === DOMnode2struct Function ===========================================
%  =======================================================================
function [s sname LeafNode] = DOMnode2struct(node, Pref, level)
[sname LeafNode] = NodeName(node);
s = [];

%% === read in node data =================================================
if (LeafNode)
  if (LeafNode>1 && ~Pref.ReadSpec), LeafNode=-1; end % tags only so ignore special nodes
  if (LeafNode>0) % supported leaf node types
    s = strtrim(char(node.getData));
    if (LeafNode==1 && Pref.Str2Num), s=str2var(s); end
  end
  if (LeafNode==3) % ProcessingInstructions need special treatment
    target = strtrim(char(node.getTarget));
    s = [target, ' ', s];
  end
  return
end
if (level>Pref.NumLevels+1), return; end

%% === read in children nodes ============================================
if (node.hasChildNodes)        % children present
  Child  = node.getChildNodes; % create array of children nodes
  nChild = Child.getLength;    % number of children

  % --- pass 1: how many children with each name -----------------------
  f = [];
  for iChild = 1:nChild        % read in each child
    [cname cLeaf] = NodeName(Child.item(iChild-1));
    if (cLeaf<0), continue; end % unsupported leaf node types
    if (~isfield(f,cname)),
      f.(cname)=0;           % initialize first time I see this name
    end
    f.(cname) = f.(cname)+1; % add to the counter
  end                        % end for iChild
  % text_nodes become CONTENT & for some reason current xmlread 'creates' a
  % lot of empty text fields so f.CONTENT value should not be trusted
  if (isfield(f,'CONTENT') && f.CONTENT>2), f.CONTENT=2; end

  % --- pass 2: store all the children ---------------------------------
  for iChild = 1:nChild        % read in each child
    [c cname cLeaf] = DOMnode2struct(Child.item(iChild-1), Pref, level+1);
    if (cLeaf && isempty(c))   % if empty leaf node than skip
      continue;                % usually empty text node or one of unhandled node types
    elseif (nChild==1 && cLeaf==1)
      s=c;                     % shortcut for a common case
    else                       % if normal node
      if (level>Pref.NumLevels), continue; end
      n = f.(cname);           % how many of them in the array so far?
      if (~isfield(s,cname))   % encountered this name for the first time
        if (n==1)              % if there will be only one of them ...
          s.(cname) = c;       % than save it in format it came in
        else                   % if there will be many of them ...
          s.(cname) = cell(1,n);
          s.(cname){1} = c;    % than save as cell array
        end
        f.(cname) = 1;         % reset the counter
      else                     % already have seen this name
        s.(cname){n+1} = c;    % add to the array
        f.(cname) = n+1;       % add to the array counter
      end
    end
  end   % for iChild
end % end if (node.hasChildNodes)

%% === Post-processing of struct's =======================================
if (isstruct(s))
  fields = fieldnames(s);
  nField = length(fields);

  % --- Post-processing: convert 'struct of arrays' to 'array of struct'
  vec = zeros(size(fields));
  for i=1:nField, vec(i) = f.(fields{i}); end
  if (numel(vec)>1 && vec(1)>1 && var(vec)==0)    % convert from struct of
    s = cell2struct(struct2cell(s), fields, 1); % arrays to array of struct
  end % if anyone knows better way to do above conversion please let me know.

  % --- Post-processing: remove special 'item' tags ---------------------
  if (isfield(s,Pref.ItemName))
    if (nField==1)
      s = s.(Pref.ItemName);         % only child: remove a level
    else
      s.CONTENT = s.(Pref.ItemName); % other children/attributes present use CONTENT
      s = rmfield(s,Pref.ItemName);
    end
  end

  % --- Post-processing: clean up CONTENT tags ---------------------
  if (isfield(s,'CONTENT'))
    if (iscell(s.CONTENT)) % && all(cellfun('isempty', s.CONTENT(2:end))))
      %msk = ~cellfun('isempty', s.CONTENT)
      %s.CONTENT = s.CONTENT(msk); % delete empty cells
      x = s.CONTENT;
      for i=length(x):-1:1, if ~isempty(x{i}), break; end; end
      if (i==1)
        s.CONTENT = x{1};   % delete cell structure
      else
        s.CONTENT = x(1:i); % delete empty cells
      end
    end
    if (nField==1)
      s = s.CONTENT;      % only child: remove a level
    end
  end
end

%% === Read in attributes ===============================================
if (node.hasAttributes && Pref.ReadAttr)
  if (~isstruct(s)),               % make into struct if is not already
    ss.CONTENT=s;
    s=ss;
  end
  Attr  = node.getAttributes;     % list of all attributes
  for iAttr = 1:Attr.getLength    % for each attribute
    name  = char(Attr.item(iAttr-1).getName);  % attribute name
    name  = str2varName(name);    % fix name if needed
    value = char(Attr.item(iAttr-1).getValue); % attribute value
    if (Pref.Str2Num), value = str2var(value); end % convert to number if possible
    s.ATTRIBUTE.(name) = value;   % save again
  end                             % end iAttr loop
end % done with attributes

%% === Post-processing: convert 'cells of structs' to 'arrays of structs'
if (isstruct(s))
  fields = fieldnames(s);     % get field names
  for iItem=1:length(s)       % for each struct in the array - usually one
    for iField=1:length(fields)
      field = fields{iField}; % get field name
      x = s(iItem).(field);
      if (iscell(x) && all(cellfun(@isstruct,x))) % it's cells of structs
        try                           % this operation fails sometimes
          s(iItem).(field) = [x{:}];  % converted to arrays of structs
        catch
          if (Pref.NoCells)
            s(iItem).(field) = forceCell2Struct(x);
          end
        end % end catch
      end
    end
  end
end

%% =======================================================================
%  === forceCell2Struct Function =========================================
%  =======================================================================
function s = forceCell2Struct(x)
% Convert cell array of structs, where not all of structs have the same
% fields, to a single array of structs

%% Convert 1D cell array of structs to 2D cell array, where each row 
% represents item in original array and each column corresponds to a unique
% field name. Array "AllFields" store fieldnames for each column
AllFields = fieldnames(x{1});     % get field names of the first struct
CellMat = cell(length(x), length(AllFields));
for iItem=1:length(x)
  fields = fieldnames(x{iItem});  % get field names of the next struct 
  for iField=1:length(fields)     % inspect all fieldnames and find those 
    field = fields{iField};       % get field name
    col = find(strcmp(field,AllFields),1);
    if isempty(col)               % no column for such fieldname yet
      AllFields = [AllFields; field];
      col = length(AllFields);    % create a new column for it
    end
    CellMat{iItem,col} = x{iItem}.(field); % store rearanged data
  end
end
%% Convert 2D cell array to array of structs
s = cell2struct(CellMat, AllFields, 2);

%% =======================================================================
%  === str2var Function ==================================================
%  =======================================================================
function val=str2var(str)
% Can this string be converted to a number? if so than do it.
val = str;
if (numel(str)==0), return; end
digits = '[Inf,NaN,pi,\t,\n,\d,\+,\-,\*,\.,e,i, ,E,I,\[,\],\;,\,]';
s = regexprep(str, digits, ''); % remove all the digits and other allowed characters
if (~all(~isempty(s)))          % if nothing left than this is probably a number
  str(strcmp(str,'\n')) = ';';  % parse data tables into 2D arrays, if any
%   try                           % try to convert to a date, like 2007-12-05
%     datenum(str);               % if successful than leave it alone
%   catch                         % if this is not a date than ...
    num = str2num(str);         % ... try converting to a number
    if(isnumeric(num) && numel(num)>0), val=num; end % if a number than save
%   end
end

%% =======================================================================
%  === str2varName Function ==============================================
%  =======================================================================
function str = str2varName(str)
% convert a sting to a valid matlab variable name
str = regexprep(str,':','_COLON_', 'once', 'ignorecase');
str = regexprep(str,'-','_DASH_'  ,'once', 'ignorecase');
if (~isvarname(str))
  str = genvarname(str);
end

%% =======================================================================
%  === NodeName Function =================================================
%  =======================================================================
function [Name LeafNode] = NodeName(node)
% get node name and make sure it is a valid variable name in Matlab.
% also get node type:
%   LeafNode=0 - normal element node,
%   LeafNode=1 - text node
%   LeafNode=2 - supported non-text leaf node,
%   LeafNode=3 - supported processing instructions leaf node,
%   LeafNode=-1 - unsupported non-text leaf node
switch (node.getNodeType)
  case node.ELEMENT_NODE
    Name = char(node.getNodeName);% capture name of the node
    Name = str2varName(Name);     % if Name is not a good variable name - fix it
    LeafNode = 0;
  case node.TEXT_NODE
    Name = 'CONTENT';
    LeafNode = 1;
  case node.COMMENT_NODE
    Name = 'COMMENT';
    LeafNode = 2;
  case node.CDATA_SECTION_NODE
    Name = 'CDATA_SECTION';
    LeafNode = 2;
  case node.DOCUMENT_TYPE_NODE
    Name = 'DOCUMENT_TYPE';
    LeafNode = 2;
  case node.PROCESSING_INSTRUCTION_NODE
    Name = 'PROCESSING_INSTRUCTION';
    LeafNode = 3;
  otherwise
    NodeType = {'ELEMENT','ATTRIBUTE','TEXT','CDATA_SECTION', ...
      'ENTITY_REFERENCE', 'ENTITY', 'PROCESSING_INSTRUCTION', 'COMMENT',...
      'DOCUMENT', 'DOCUMENT_TYPE', 'DOCUMENT_FRAGMENT', 'NOTATION'};
    Name = char(node.getNodeName);% capture name of the node
    warning('xml_io_tools:read:unkNode', ...
      'Unknown node type encountered: %s_NODE (%s)', NodeType{node.getNodeType}, Name);
    LeafNode = -1;
end

%% =======================================================================
%  === mergeStructs Function =================================================
%  =======================================================================
function res = mergeStructs(x,y)
% From: http://stackoverflow.com/a/6271161
if isstruct(x) && isstruct(y)
    res = x;
    names = fieldnames(y);
    for fnum = 1:numel(names)
        if isfield(x,names{fnum})
            res.(names{fnum}) = mergeStructs(x.(names{fnum}),y.(names{fnum}));
        else
            res.(names{fnum}) = y.(names{fnum});
        end
    end
else
    res = y;
end