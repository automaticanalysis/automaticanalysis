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

switch nargin
    case {0 1}
        aa_info = aaClass('nopath', 'nogreet');
        defaultparameters = aa_info.parameter_xml_filename;
        if nargin
            tasklistxml = varargin{1};
        else
            warning('You must provide a tasklist to aarecipe.\n');
        end
    case 2
        defaultparameters = varargin{1};
        tasklistxml = varargin{2};
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
    assert(~strcmp(resp,'No (Exit)'), 'exiting'); % TODO: This does not take into account the cancel button / dismissing the dialog / random userinput. Compare to the positive option only!
    defaultparameters = aas_create_parameter_xml(defaultparameters, isGUI);
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

%% pruneEmptyBranches
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

%% checkhasrequiredfields
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

%% aa_unique
function [label, index] = aa_unique(vals)

try
    [label, index] = unique(vals, 'legacy');
catch
    [label, index] = unique(vals);
end

end
