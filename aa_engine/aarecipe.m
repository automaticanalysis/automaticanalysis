% Load parameter defaults and tasklist into the structure "aap"
%
% FORMAT aap = aarecipe(tasklist)
% Parameter defaults are loaded from ~/.aa/aap_parameters_user.xml
%   - tasklist: XML-file containing the list of modules
%
% FORMAT aap = aarecipe(parameters,tasklist)
%   - parameters: XML-file containing the parameter defaults
%   - tasklist: XML-file containing the list of modules
%
% If the specified parameters file does not exist, the user is prompted to generate one
% by copying one of the files that are included in AA. If you use this functionality you
% may need to make further manual edits to the parameter file.
%
% Rhodri Cusack
% Tibor Auer MRC CBU Cambridge 2016

function aap = aarecipe(varargin)

isGUI = ~any(strcmp(varargin,'nogui'));
varargin(strcmp(varargin,'nogui')) = [];

Pref.ReadAttr=0;

switch(numel(varargin))
    case 0
        warning('You must provide a tasklist to aarecipe.\n');
    case 1
        defaultparameters='aap_parameters_user.xml';
        tasklistxml=varargin{1};
    case 2
        defaultparameters=varargin{1};
        tasklistxml=varargin{2};
end

clear aap

% First work on default parameters
if ~exist(defaultparameters,'file')
    fprintf(...
        'Cannot find parameters file %s, opening user interface to generate a new file\n',...
        defaultparameters);
    resp = userinput('questdlg',...
        sprintf('Cannot find parameters file %s\nSeed new parameter file from existing default?',defaultparameters), ...
        'Parameter file', ...
        'Yes','No (Exit)','No (Exit)','GUI',isGUI);
    assert(~strcmp(resp,'No (Exit)'), 'exiting');
    % if we made it here, we are seeding a new parameters file
    % we have default parameters
    defaultdir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'aa_parametersets');
    [seedparam, rootpath] = userinput('uigetfile',{'*.xml','All Parameter Files' },'Desired seed parameter',defaultdir,'GUI',isGUI);
    assert(ischar(seedparam), 'exiting');
    seedparam = fullfile(rootpath, seedparam);

    % initialise the save dialogue in the current aap.acq_details.root if specified
    xml=xml_read(seedparam);
    configdir = fullfile(getenv('HOME'),'.aa');
    % generate new parameters file N.B.: in networks with shared resources
    % average user may not be able to write into aa_paremetersets
    [defaultparameters, rootpath] = userinput('uiputfile',{'*.xml','All Parameter Files' },...
        'Location of the parameters file',fullfile(configdir, defaultparameters),'GUI',isGUI);
    assert(ischar(defaultparameters), 'exiting');
    destination = fullfile(rootpath, defaultparameters);

    analysisroot = aas_expandpathbyvars(xml.acq_details.root.CONTENT);
    aas_makedir([], analysisroot);
    analysisroot = userinput('uigetdir',analysisroot,'Location of analyses by default','GUI',isGUI);

    create_minimalXML(seedparam, destination, analysisroot);
    assert(exist(destination,'file')>0,'failed to create %s', defaultparameters);

    % N.B. we don't actually modify defaultparameters - it should now be on the path. But
    % let's double check. It might not be e.g. if you haven't actually added AA to your
    % path properly before calling this function.
    assert(exist(defaultparameters,'file')>0, ...
        'could not find %s - Are you sure it is on your path?',...
        defaultparameters);

    if isGUI
        h = msgbox(sprintf('New parameter set in %s has been created.\nYou may need to edit this file further to reflect local configuration.',destination),'New parameters file','Warn');
        waitfor(h);
    else
        fprintf('\nNew parameter set in %s has been created.\nYou may need to edit this file further to reflect local configuration.\n',destination);
    end
end

schema = xml_read(defaultparameters);
aap = processattributes(schema);
aap.schema = schema;

% And now load up task list
if exist('tasklistxml','var')
    if ~exist(tasklistxml,'file')
        aas_log(aap,true,sprintf('Cannot find file %s as specified in call to aarecipe',tasklistxml));
    end
    Pref.ReadAttr=1;
    Pref.ReadSpec=0;
    Pref.ReadAttr=0;
    xml_tasklist=xml_read(tasklistxml,Pref);
    aap=struct_update(aap,xml_tasklist);
    % Now load up each modules parameters
    aap.tasksettings=[];


    % Now load up parameters for each of the modules to be used
    % For initialisation modules...
    for task=1:length(aap.tasklist.initialisation.module)
        [aap, index]=aas_addtaskparameters(aap,aap.tasklist.initialisation.module(task).name);
        aap.tasklist.initialisation.module(task).index=index;
    end
    % ...and for main modules...


    % Deal with branches - only supported for main modules
    aap.tasklist.main=processbranch(aap,aap.directory_conventions.analysisid_suffix,'*',aap.tasklist.main,1);

    aap.tasklist.main.module = pruneEmptyBranches(aap.tasklist.main.module);

    % These fields were only used for creating the branching structures
    aap.tasklist.main.module = rmfield(aap.tasklist.main.module, {'branchID', 'ignorebranches'});

    % Calculate task parameters and indices for modules
    % e.g., aamod_smooth_01, aamod_smooth_02 etc
    if isfield(aap.tasklist.main,'module')
        for task=1:length(aap.tasklist.main.module)
            [aap, index]=aas_addtaskparameters(aap,aap.tasklist.main.module(task).name,aap.tasklist.main.module(task).aliasfor);
            aap.tasklist.main.module(task).index=index;
            if isfield(aap.schema.tasksettings.(aap.tasklist.main.module(task).name)(index).ATTRIBUTE,'mfile_alias')
                aap.tasklist.main.module(task).aliasfor = aap.schema.tasksettings.(aap.tasklist.main.module(task).name)(index).ATTRIBUTE.mfile_alias;
            end
        end
    end

    % When processing the branches, the indices for modules repeated had
    % not yet been calculated
    % Now we've got the stage tags, we can go through and adorn the
    % previously numeric dependencies
    for task=1:length(aap.tasklist.main.module)
        tbcf_num=aap.tasklist.main.module(task).tobecompletedfirst;
        tbcf_cell={};
        for tbcfind=1:length(tbcf_num)
            if (isnumeric(tbcf_num(tbcfind)))
                tbcf_cell{tbcfind}=aas_getstagetag(aap,tbcf_num(tbcfind));
            end
        end
        aap.tasklist.main.module(task).tobecompletedfirst=tbcf_cell;
    end

end

% Make copy of aap
aap.aap_beforeuserchanges=[];
aap.aap_beforeuserchanges=aap;

end

%% Recursively process nodes
function node = processattributes(node)
if isstruct(node)
    if isfield(node,'ATTRIBUTE')
        if isfield(node,'CONTENT')
            attr = node.ATTRIBUTE;
            node = node.CONTENT;
            if isfield(attr,'ui')
                switch attr.ui
                    case {'text' 'dir' 'dir_allowwildcards' 'dir_part_allowwildcards' 'dir_part' 'file'}
                        node = char(node);
                        % TODO
                        %                     case {'dir_list','optionlist'}
                        %                     case {'structarray'}
                        %                     case {'intarray' 'rgb'}
                        %                     case {'double'}
                        %                     case {'int'}
                        %                     case {'yesno'}
                end
            end
            return
        else
            node = rmfield(node,'ATTRIBUTE');
        end
    end
    for f = fieldnames(node)'
        node.(f{1}) = cell2mat(arrayfun(@(x) processattributes(x), node.(f{1}), 'UniformOutput', false)); % deal with arrays
    end
end
end

%% Recursively process branches
% Used even for the "top level" branch always present in every
% tasklist
function [outstages]=processbranch(aap,analysisid_suffix,selected_sessions,branch,branchID)

if (~isfield(branch,'module') || isempty(branch.module))
    %     aas_log(aap,false,sprintf('The branch in your tasklist with analysis id suffix %s appears to be empty. Have you remembered the <module> tag that must surround the <name> tag of each module, e.g., <module><name>aamod_smooth</name></module>?',analysisid_suffix));
    %     outstages=[];
    extrastages.module.name = 'emptybranch';
    extrastages.module.extraparameters.aap.directory_conventions.analysisid_suffix=analysisid_suffix;
    extrastages.module.extraparameters.aap.acq_details.selected_sessions=selected_sessions;
    extrastages.module.branchID = branchID;
    extrastages.module.tobecompletedfirst = 0;
    extrastages = checkhasrequiredfields(extrastages);
    outstages = extrastages;

else

    % For each module
    for stagenum=1:length(branch.module)
        % Is it a regular module...
        if (~isfield(branch.module(stagenum),'branch') || isempty(branch.module(stagenum).branch))
            extrastages=[];
            extrastages.module=branch.module(stagenum);
            extrastages.module.extraparameters.aap.directory_conventions.analysisid_suffix=analysisid_suffix;
            extrastages.module.extraparameters.aap.acq_details.selected_sessions=selected_sessions;
            extrastages.module.tobecompletedfirst = 0;
            extrastages.module.branchID = branchID;

            %...or a set of branches
        else
            clear extrastages
            selected_sessions_main = selected_sessions; % backup selected_sessions
            for branchnum=1:length(branch.module(stagenum).branch)
                try
                    analysisid_suffix_append=branch.module(stagenum).branch(branchnum).analysisid_suffix;
                catch
                    analysisid_suffix_append='';
                end
                try
                    selected_sessions=strtrim(branch.module(stagenum).branch(branchnum).selected_sessions);
                    if (isempty(selected_sessions))
                        selected_sessions='*';
                    end
                catch exception
                    selected_sessions='*';
                end

                [extrastages_added] = processbranch(aap, [analysisid_suffix analysisid_suffix_append],selected_sessions,branch.module(stagenum).branch(branchnum),1);


                if (~isempty(extrastages_added))

                    % And add the new stages to our list
                    if (exist('extrastages','var'))

                        numPrevBranchStages = length(extrastages.module);

                        dependI = find([extrastages_added.module.tobecompletedfirst]~=0);   % Stages that have dependencies within this set of new stages
                        noDependI = find([extrastages_added.module.tobecompletedfirst]==0); % Stages that don't have dependencies

                        % Update pointers to account for stages in the previous in branch
                        if ~isempty(dependI)
                            extrastages_added.module(dependI) = arrayfun(@(x) setfield(x, 'tobecompletedfirst', x.tobecompletedfirst + numPrevBranchStages), extrastages_added.module(dependI));
                        end

                        if (~isempty(noDependI))
                            extrastages_added.module(noDependI) = arrayfun(@(x) setfield(x, 'tobecompletedfirst', 0), extrastages_added.module(noDependI));
                        end

                        % Update the branchIDs in the current branch to account for the number in the previous branch number of branches in the previous branch
                        numPrevBranches = length(aa_unique([extrastages.module.branchID]));
                        extrastages_added.module = arrayfun(@(x) setfield(x, 'branchID', x.branchID + numPrevBranches), extrastages_added.module);

                        extrastages.module=[extrastages.module extrastages_added.module];
                    else
                        extrastages.module=extrastages_added.module;
                    end

                end % End if (~isempty(extrastages_added))

            end % End for branchnum=1:length(....

            selected_sessions = selected_sessions_main; % restore selected_sessions
        end

        if (isfield(extrastages.module,'branch'))
            extrastages.module = rmfield(extrastages.module,'branch');
        end

        extrastages = checkhasrequiredfields(extrastages);
        if (exist('outstages','var'))

            % Get a list of all the branch IDs that exist in the output
            oBranchIDs = [outstages.module.branchID];
            [oLabels, oIndex] = aa_unique(oBranchIDs);
            numOutBranches = length(oLabels);

            % Get a list of all the branch IDs that exist in the new stuff
            nBranchIDs = [extrastages.module.branchID];
            [nLabels, ~] = aa_unique(nBranchIDs);
            numNewBranches = length(nLabels);

            % Number of stages in the new stuff
            numNewStages = length(extrastages.module);

            % Basically, we repeat the extrastages and them to each branch that exists in the output.
            for oB = 1 : numOutBranches
                newStages = extrastages;

                % Update the dependencies: some stages will have no
                % dependencies (e.g., the first stage of a new branch),
                % they will point to the last stage of the output branch.
                noDependI = find([newStages.module.tobecompletedfirst]==0);
                dependI = find([newStages.module.tobecompletedfirst]~=0);
                if ~isempty(noDependI)
                    newStages.module(noDependI) = arrayfun(@(x) setfield(x, 'tobecompletedfirst', oIndex(oB)), newStages.module(noDependI));
                end

                % Or, some new stages are dependent on other new stages, so
                % they will have to be updated to account for the repeating
                % nature of this operation.
                if ~isempty(dependI)
                    newStages.module(dependI) = arrayfun(@(x) setfield(x, 'tobecompletedfirst', x.tobecompletedfirst+length(outstages.module)), newStages.module(dependI));
                end

                % Update the analysis suffixes
                outBranchSuffix = outstages.module(oIndex(oB)).extraparameters.aap.directory_conventions.analysisid_suffix;
                if ~isempty(outBranchSuffix)
                    for stage = 1 : numNewStages
                        newStagesSuffix = newStages.module(stage).extraparameters.aap.directory_conventions.analysisid_suffix;
                        newStages.module(stage).extraparameters.aap.directory_conventions.analysisid_suffix = [outBranchSuffix newStagesSuffix(length(analysisid_suffix)+1:end)];
                    end
                end

                % Update the selected_sessions
                newStagestoRemove = [];
                for stage = 1 : numNewStages
                    outBranchSel = outstages.module(oIndex(oB)).extraparameters.aap.acq_details.selected_sessions;
                    newStagesSel = newStages.module(stage).extraparameters.aap.acq_details.selected_sessions;

                    outBranchSelC = textscan(outBranchSel,'%s','delimiter',' '); outBranchSelC = outBranchSelC{1};
                    newStagesSelC = textscan(newStagesSel,'%s','delimiter',' '); newStagesSelC = newStagesSelC{1};
                    sessionsSelC = textscan(selected_sessions,'%s','delimiter',' '); sessionsSelC = sessionsSelC{1};
                    if any(strcmp(outBranchSelC,'*'))
                        outBranchSel = union(newStagesSelC,sessionsSelC);
                    else
                        outBranchSel = outBranchSelC;
                    end
                    if any(strcmp(newStagesSelC,'*'))
                        newStagesSel = union(outBranchSelC,sessionsSelC);
                    else
                        newStagesSel = newStagesSelC;
                    end
                    if any(strcmp(sessionsSelC,'*'))
                        sessionsSel = union(newStagesSelC,outBranchSelC);
                    else
                        sessionsSel = sessionsSelC;
                    end

                    sessSel = intersect(outBranchSel,intersect(newStagesSel,sessionsSel));
                    if isempty(sessSel)
                        newStagestoRemove(end+1) = stage;
                        continue;
                    end
                    sessSel = sprintf('%s ',sessSel{:}); sessSel(end) = '';

                    newStages.module(stage).extraparameters.aap.acq_details.selected_sessions = sessSel;
                end
                newStages.module(newStagestoRemove) = [];

                % Update the branchIDs
                newStages.module = arrayfun(@(x) setfield(x, 'branchID', x.branchID+(oB-1)*numNewBranches), newStages.module);

                % Add stages
                outstages.module = [outstages.module newStages.module];

            end % End for oB = 1 : numOutputBranches

        else
            outstages.module = extrastages.module;
        end

    end % End for loop stagenum

end
end % End Function processbranch

function outstages = pruneEmptyBranches(outstages)

% Indices of empty branches
eBranchInd = find(strcmp({outstages.name}, 'emptybranch'));

for m = 1 : length(outstages)
    if isstr(outstages(m).ignorebranches)
        bNames = regexp(outstages(m).ignorebranches, '\w+', 'match');
        if find(strcmpi(bNames, outstages(m).extraparameters.aap.directory_conventions.analysisid_suffix))
            eBranchInd(end+1) = m;
        end
    end
end

eBranchInd = sort(eBranchInd);

% Redirect all modules that point to an empty branch to instead point
% to the module that the empty branch points to.
for e = 1 : length(eBranchInd)
    mInd = find([outstages.tobecompletedfirst] == eBranchInd(e));
    outstages(mInd) = arrayfun(@(x) setfield(x, 'tobecompletedfirst', outstages(eBranchInd(e)).tobecompletedfirst), outstages(mInd));
end

% Now remove the 'emptybranch' modules, and update all the pointers to
% account for the missing module. Working backwards ensures that we don't
% mess up if we shift, then remove, then shift again.
for e = length(eBranchInd) : -1 : 1
    mInd = find([outstages.tobecompletedfirst] > eBranchInd(e));
    outstages(mInd) = arrayfun(@(x) setfield(x, 'tobecompletedfirst', x.tobecompletedfirst-1), outstages(mInd));
    outstages(eBranchInd(e)) = [];
end

% Get rid of any pointers that are 0
mInd = find([outstages.tobecompletedfirst] == 0);
outstages(mInd) = arrayfun(@(x) setfield(x, 'tobecompletedfirst', []), outstages(mInd));

end

function outstages=checkhasrequiredfields(outstages)

% Check all required fields are present
reqflds={'index','name','epiprefix','tobecompletedfirst','extraparameters','aliasfor','remotestream','branchID', 'ignorebranches'};

for stagenum=1:length(outstages.module)
    for reqfldsind=1:length(reqflds)
        if (~isfield(outstages.module(stagenum),reqflds{reqfldsind}))
            outstages.module(stagenum).(reqflds{reqfldsind})=[];
        end
    end
end

end

function [label, index] = aa_unique(vals)

try
    [label, index] = unique(vals, 'legacy');
catch
    [label, index] = unique(vals);
end

end

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

function varargout = userinput(varargin)
% Examples:
% resp = userinput('questdlg',sprintf('Cannot find parameters file %s\nSeed new parameter file from existing default?','paramfile'), 'Parameter file', 'Yes','No (Exit)','No (Exit)','GUI',true);
% [seedparam, rootpath] = userinput('uigetfile',{'*.xml','All Parameter Files' },'Desired seed parameter',defaultdir,'GUI',true);
% [defaultparameters, rootpath] = userinput('uiputfile',{'*.xml','All Parameter Files' }, 'Location of the parameters file and analyses by default',fullfile(pwd,defaultparameters),'GUI',true);

isGUI = true;
iParam = find(strcmpi(varargin,'gui'),1);
if ~isempty(iParam)
    isGUI = varargin{iParam+1};
    varargin(iParam:iParam+1) = [];
end

switch varargin{1}
    case 'questdlg' % question, title, btn1, btn2, default
        if isGUI
            varargout{1} = questdlg(varargin{2:end});
        else
            btns = varargin(4:end-1);
            msgBtn = sprintf(' / %s',btns{:}); msgBtn(1:3) = '';
            respList = cellfun(@(x) lower(strtok(x)), btns, 'UniformOutput', false);
            while true
                resp = input([varargin{2} ' (' msgBtn '):' ],'s');
                resp = btns(cellfun(@(x) strcmp(resp,x) || (resp==x(1)), respList));
                if ~isempty(resp), break; end
            end
            varargout{1} = resp{1};
        end
    case  'uigetdir'
        if isGUI
            varargout{1} = uigetdir(varargin{2:end});
        else
            rootpath = input([varargin{3} ' (or leave empty for ' varargin{2} '):'],'s');
            if isempty(rootpath), rootpath = varargin{2}; end

            varargout{1} = rootpath;
        end
    case  'uigetfile'
        if isGUI
            [varargout{1}, varargout{2}] = uigetfile(varargin{2:end});
        else
            defaultdir = varargin{4};
            defaultnames = dir(fullfile(defaultdir,varargin{2}{1}));
            fprintf('%s in %s:\n', varargin{2}{2}, defaultdir);
            fprintf('%s\n',defaultnames.name);
            while true
                seedparam = input([varargin{3} ' (or leave empty to abort):'],'s');
                % filter out extension (so we are robust to whether this is provided or not)
                [rootpath,seedparam,~] = fileparts(seedparam);
                if isempty(rootpath), rootpath = defaultdir; end
                seedparam = [seedparam varargin{2}{1}(2:end)];

                if exist(fullfile(rootpath,seedparam),'file'), break;
                else
                    fprintf('Could not find file %s. Please try again!\n',seedparam);
                end
            end

            varargout{1} = seedparam;
            varargout{2} = rootpath;
        end
    case  'uiputfile'
        if isGUI
            [varargout{1}, varargout{2}] = uiputfile(varargin{2:end});
        else
            [defaultdir, defaultseed]= fileparts(varargin{4});
            seedparam = input([varargin{3} ' (or leave empty for ' varargin{4} '):'],'s');
            % filter out extension (so we are robust to whether this is provided or not)
            [rootpath,seedparam] = fileparts(seedparam);
            if isempty(rootpath), rootpath = defaultdir; end
            if isempty(seedparam), seedparam = defaultseed; end

            varargout{1} = [seedparam varargin{2}{1}(2:end)];
            varargout{2} = rootpath;
        end
    otherwise
        error('Function %s is not an existing function or not implemented!',varargin{1});
end
end
