function [aap]=aarecipe(varargin)

switch(nargin)
    case 0
        fprintf('You must provide a tasklist to aarecipe.\n');
    case 1
        defaultparameters='aap_parameters_defaults.xml';
        tasklistxml=varargin{1};
    case 2
        defaultparameters=varargin{1};
        tasklistxml=varargin{2};
end;

clear aap

% First work on default parameters
if (~length(which(defaultparameters)))
    fprintf('Cannot find file %s as specified in call to aarecipe\n',defaultparameters);
end;

Pref.ReadAttr=0;
aap=xml_read(defaultparameters,Pref);
aap.schema=xml_read(defaultparameters);


% And now load up task list
if (exist('tasklistxml','var'))
    if (~length(which(tasklistxml)))
        aas_log(aap,true,sprintf('Cannot find file %s as specified in call to aarecipe',tasklistxml));
    end;
    Pref.ReadAttr=1;
    Pref.ReadSpec=0;
    xml_tasklist.schema=xml_read(tasklistxml,Pref);
    Pref.ReadAttr=0;
    xml_tasklist=xml_read(tasklistxml,Pref);
    aap=catstruct(aap,xml_tasklist);
    % Now load up each modules parameters
    aap.tasksettings=[];
    
    
    % Now load up parameters for each of the modules to be used
    % For initialisation modules...
    for task=1:length(aap.tasklist.initialisation.module)
        [aap index]=aas_addtaskparameters(aap,aap.tasklist.initialisation.module(task).name);
        aap.tasklist.initialisation.module(task).index=index;
    end;
    % ...and for main modules...
    
    
    % Deal with branches - only supported for main modules
    aap.tasklist.main=processbranch(aap,aap.directory_conventions.analysisid_suffix,'*',aap.tasklist.main,1);
    
    aap.tasklist.main.module = pruneEmptyBranches(aap.tasklist.main.module);
    
<<<<<<< HEAD
    aap.tasklist.main.module = rmfield(aap.tasklist.main.module, {'branchID', 'branchDepth', 'ignorebranches'});
=======
    aap.tasklist.main.module = rmfield(aap.tasklist.main.module, {'branchID', 'ignorebranches'});
>>>>>>> c71b509be7b591da4669e7e639ba43c795cc5ebf
    
    % Calculate task parameters and indices for modules
    % e.g., aamod_smooth_01, aamod_smooth_02 etc
    if isfield(aap.tasklist.main,'module')
        for task=1:length(aap.tasklist.main.module)
            [aap index]=aas_addtaskparameters(aap,aap.tasklist.main.module(task).name,aap.tasklist.main.module(task).aliasfor);
            aap.tasklist.main.module(task).index=index;
        end;
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
            end;
        end;
        aap.tasklist.main.module(task).tobecompletedfirst=tbcf_cell;
    end;
    
end;

% And copy in SPM defaults
global defaults
try
    aap.spm.defaults=spm_get_defaults;
catch
    
    if (~isstruct(defaults))
        aas_log(aap,true,'You must launch SPM before running aa- try typing spm fmri;');
    else
        aap.spm.defaults=defaults;
    end;
end;

% Make copy of aap
aap.aap_beforeuserchanges=[];
aap.aap_beforeuserchanges=aap;

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
<<<<<<< HEAD
    extrastages.module.branchDepth = 1; %branchDepth;
    
=======
    extrastages.module.tobecompletedfirst = 0;
>>>>>>> c71b509be7b591da4669e7e639ba43c795cc5ebf
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
<<<<<<< HEAD
            extrastages.module.branchDepth = 1; % branchDepth;
=======
>>>>>>> c71b509be7b591da4669e7e639ba43c795cc5ebf
            
            %...or a set of branches
        else
            clear extrastages
            
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
                
<<<<<<< HEAD
                [extrastages_added] = processbranch(aap,[analysisid_suffix analysisid_suffix_append],selected_sessions,branch.module(stagenum).branch(branchnum),1);
=======
                [extrastages_added] = processbranch(aap, [analysisid_suffix analysisid_suffix_append],selected_sessions,branch.module(stagenum).branch(branchnum),1);
>>>>>>> c71b509be7b591da4669e7e639ba43c795cc5ebf
                
                
                if (~isempty(extrastages_added))
                    
                    % And add the new stages to our list
                    if (exist('extrastages','var'))
                        
                        numPrevBranchStages = length(extrastages.module);
                        
                        dependI = find([extrastages_added.module.tobecompletedfirst]~=0);   % Stages that have dependencies within this set of new stages
                        noDependI = find([extrastages_added.module.tobecompletedfirst]==0); % Stages that don't have dependencies
                        
                        % Update pointers to account for stages in the previous in branch
                        extrastages_added.module(dependI) = arrayfun(@(x) setfield(x, 'tobecompletedfirst', x.tobecompletedfirst + numPrevBranchStages), extrastages_added.module(dependI));
                        extrastages_added.module(noDependI) = arrayfun(@(x) setfield(x, 'tobecompletedfirst', 0), extrastages_added.module(noDependI));
                        
                        % Update the branchIDs in the current branch to account for the number in the previous branch number of branches in the previous branch
                        numPrevBranches = length(unique([extrastages.module.branchID]));
                        
                        extrastages_added.module = arrayfun(@(x) setfield(x, 'branchID', x.branchID + numPrevBranches), extrastages_added.module);
                        
                        extrastages.module=[extrastages.module extrastages_added.module];
                    else
                        extrastages.module=extrastages_added.module;
                    end
                    
                end % End if (~isempty(extrastages_added))
                
            end % End for branchnum=1:length(....
            
        end
        
        if (isfield(extrastages.module,'branch'))
            extrastages.module = rmfield(extrastages.module,'branch');
        end;
        
        extrastages = checkhasrequiredfields(extrastages);
        if (exist('outstages','var'))
            
            % Get a list of all the branch IDs that exist in the output
            oBranchIDs = [outstages.module.branchID];
            [oLabels oIndex] = unique(oBranchIDs);
            numOutBranches = length(oLabels);
            
            % Get a list of all the branch IDs that exist in the new stuff
            nBranchIDs = [extrastages.module.branchID];
            [nLabels nIndex] = unique(nBranchIDs);
            numNewBranches = length(nLabels);
            
            % Number of stages in the new stuff
            numNewStages = length(extrastages.module);
            
<<<<<<< HEAD
            if numInputBranches > 1
                maxPrevBranchDepth = max([outstages.module.branchDepth]);
                extrastages.module = arrayfun(@(x) setfield(x, 'branchDepth', maxPrevBranchDepth + x.branchDepth), extrastages.module);
            end
            
            % Basically, we 'repmat' the extrastages to add them to
            % each branch that exists in the output.
            branchRepeatedStages = struct([]);
=======
            % Basically, we repeat the extrastages and them to each branch that exists in the output.
>>>>>>> c71b509be7b591da4669e7e639ba43c795cc5ebf
            for oB = 1 : numOutBranches
                newStages = extrastages;
                
<<<<<<< HEAD
                outBranchSuffix = outstages.module(oIndex(oB)).extraparameters.aap.directory_conventions.analysisid_suffix;
                
                % Create a copy of all the modules to be added
                branchStages = extrastages;
                
                % For each branch in the new stages...
                for iB = 1 : numInputBranches
                    
                    % Get the indices of stages that belong to this branch
                    iBInd = find(iBranchIDs == iLabels(iB));
                    
                    % Now, for each of these stages update the IDs, path
                    % suffix, and dependency
                    for m = 1 : length(iBInd)
                        
                        % Set the branch ID
                        branchStages.module(iBInd(m)).branchID = (oB-1)*numInputBranches + iB + branchID-1;
                        
                        % Concatenate the path suffixes, but only if we are connecting up different branches
                        thisBranchSuffix = branchStages.module(iBInd(m)).extraparameters.aap.directory_conventions.analysisid_suffix;
                        %                         if ~strcmp(outBranchSuffix, thisBranchSuffix)
                        %                             branchStages.module(iBInd(m)).extraparameters.aap.directory_conventions.analysisid_suffix = strcat(outBranchSuffix, thisBranchSuffix);
                        %                         end
                        
                        if (outstages.module(oIndex(oB)).branchDepth <= branchStages.module(iBInd(m)).branchDepth) && ~strcmp(outBranchSuffix, thisBranchSuffix)
                            branchStages.module(iBInd(m)).extraparameters.aap.directory_conventions.analysisid_suffix = strcat(outBranchSuffix, thisBranchSuffix);
                        else
                            branchStages.module(iBInd(m)).extraparameters.aap.directory_conventions.analysisid_suffix = outBranchSuffix;
                        end
                        
                        % Set the dependencies... gets a little tricky
                        
                        % Connect the first stage of the new branch to the
                        % last stage of the current output branch.
                        if m == 1
                            branchStages.module(iBInd(m)).tobecompletedfirst = oIndex(oB);
                            % But we might already point to a previous
                            % stage, especicially if we are using lots of
                            % sub branches. So correct for that...
                        elseif isfield(branchStages.module(iBInd(m)), 'tobecompletedfirst') && ~isempty(branchStages.module(iBInd(m)).tobecompletedfirst)
                            branchStages.module(iBInd(m)).tobecompletedfirst = branchStages.module(iBInd(m)).tobecompletedfirst + length(branchRepeatedStages) + length(outstages.module);
                        end
                        
                    end % End for m = 1 : length(iBInd)
                    
                end % End for iB = 1 : numInputBranches
=======
                % Update the dependencies: some stages will have no
                % dependencies (e.g., the first stage of a new branch),
                % they will point to the last stage of the output branch.
                noDependI = find([newStages.module.tobecompletedfirst]==0);
                dependI = find([newStages.module.tobecompletedfirst]~=0);
                newStages.module(noDependI) = arrayfun(@(x) setfield(x, 'tobecompletedfirst', oIndex(oB)), newStages.module(noDependI));
                
                % Or, some new stages are dependent on other new stages, so
                % they will have to be updated to account for the repeating
                % nature of this operation.
                newStages.module(dependI) = arrayfun(@(x) setfield(x, 'tobecompletedfirst', x.tobecompletedfirst+length(outstages.module)), newStages.module(dependI));
                
                % Update the analysis suffixes
                outBranchSuffix = outstages.module(oIndex(oB)).extraparameters.aap.directory_conventions.analysisid_suffix;
                if ~isempty(outBranchSuffix)
                    for stage = 1 : numNewStages
                        newStagesSuffix = newStages.module(stage).extraparameters.aap.directory_conventions.analysisid_suffix;
                        newStages.module(stage).extraparameters.aap.directory_conventions.analysisid_suffix = [outBranchSuffix newStagesSuffix(length(analysisid_suffix)+1:end)];
                    end
                end
                
                % Update the branchIDs
                newStages.module = arrayfun(@(x) setfield(x, 'branchID', x.branchID+(oB-1)*numNewBranches), newStages.module);
                
                
                outstages.module = [outstages.module newStages.module];
>>>>>>> c71b509be7b591da4669e7e639ba43c795cc5ebf
                
                
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
<<<<<<< HEAD
   if isstr(outstages(m).ignorebranches)
       bNames = regexp(outstages(m).ignorebranches, '\w+', 'match');
        if find(strcmpi(bNames, outstages(m).extraparameters.aap.directory_conventions.analysisid_suffix))
            eBranchInd(end+1) = m;
        end
   end
=======
    if isstr(outstages(m).ignorebranches)
        bNames = regexp(outstages(m).ignorebranches, '\w+', 'match');
        if find(strcmpi(bNames, outstages(m).extraparameters.aap.directory_conventions.analysisid_suffix))
            eBranchInd(end+1) = m;
        end
    end
>>>>>>> c71b509be7b591da4669e7e639ba43c795cc5ebf
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
<<<<<<< HEAD
reqflds={'index','name','epiprefix','tobecompletedfirst','extraparameters','aliasfor','remotestream','branchID','branchDepth', 'ignorebranches'};
=======
reqflds={'index','name','epiprefix','tobecompletedfirst','extraparameters','aliasfor','remotestream','branchID', 'ignorebranches'};
>>>>>>> c71b509be7b591da4669e7e639ba43c795cc5ebf

for stagenum=1:length(outstages.module)
    for reqfldsind=1:length(reqflds)
        if (~isfield(outstages.module(stagenum),reqflds{reqfldsind}))
            outstages.module(stagenum).(reqflds{reqfldsind})=[];
        end;
    end;
end

end
