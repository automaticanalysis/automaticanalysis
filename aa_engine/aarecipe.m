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
    
    aap.tasklist.main.module = rmfield(aap.tasklist.main.module, {'branchID'});
    
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
    aas_log(aap,false,sprintf('The branch in your tasklist with analysis id suffix %s appears to be empty. Have you remembered the <module> tag that must surround the <name> tag of each module, e.g., <module><name>aamod_smooth</name></module>?',analysisid_suffix));
    outstages=[];
else

    % For each module
    for stagenum=1:length(branch.module)
        % Is it a regular module...
        if (~isfield(branch.module(stagenum),'branch') || isempty(branch.module(stagenum).branch))
            extrastages=[];
            extrastages.module=branch.module(stagenum);
            extrastages.module.extraparameters.aap.directory_conventions.analysisid_suffix=analysisid_suffix;
            extrastages.module.extraparameters.aap.acq_details.selected_sessions=selected_sessions;
            extrastages.module.branchID = branchID;
            
            %...or a set of branches
        else
            clear extrastages
%             if (stagenum>1)
%                 branchdependency=stagenum-1;
%             else
%                 branchdependency='[none]';
%             end;
            
            numSubBranches = 0;
            for branchnum=1:length(branch.module(stagenum).branch)
                try
                    analysisid_suffix_append=branch.module(stagenum).branch(branchnum).analysisid_suffix;
                catch
                    analysisid_suffix_append='';
                end;
                try
                    selected_sessions=strtrim(branch.module(stagenum).branch(branchnum).selected_sessions);
                    if (isempty(selected_sessions))
                        selected_sessions='*';
                    end;
                catch exception
                    selected_sessions='*';
                end;
                
                numSubBranches = numSubBranches + 1;
                [extrastages_added] = processbranch(aap,[analysisid_suffix analysisid_suffix_append],selected_sessions,branch.module(stagenum).branch(branchnum),numSubBranches);              
                
                
                if (~isempty(extrastages_added))
%                     % Set dependency of first stage of branch to last stage of
%                     % prededing analysis
%                     if (isempty(extrastages_added.module(1).tobecompletedfirst))
%                         extrastages_added.module(1).tobecompletedfirst=branchdependency;
%                     end;
%                     
                    % And add the new stages to our list
                    if (exist('extrastages','var'))
                        
                        for m = 1 : length(extrastages_added.module)
                            extrastages_added.module(m).tobecompletedfirst = extrastages_added.module(m).tobecompletedfirst + length(extrastages.module);
                        end
                            
                        extrastages.module=[extrastages.module extrastages_added.module];
                    else
                        extrastages.module=extrastages_added.module;
                    end;
                end;
            end;
            
        end;
        
        
        if (isfield(extrastages.module,'branch'))
            extrastages.module=rmfield(extrastages.module,'branch');
        end;
        extrastages=checkhasrequiredfields(extrastages);
        if (exist('outstages','var'))
            
            % Get a list of all the branch IDs that exist in the output
            oBranchIDs = [outstages.module.branchID];
            [oLabels oIndex] = unique(oBranchIDs);
            numOutBranches = length(oLabels);
            
            % Get a list of all the branch IDs in the new stages
            iBranchIDs = [extrastages.module.branchID];
            [iLabels iIndex] = unique(iBranchIDs);
            numInputBranches = length(iLabels);           
            
            % Basically, we 'repmat' the extrastages to add them to
            % each branch that exists in the output.
            branchRepeatedStages = struct([]);
            for oB = 1 : numOutBranches
                
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
                        
                        % Concatenate the path suffixes, but only if we are
                        % connecting up different branches
                        outBranchSuffix = outstages.module(oIndex(oB)).extraparameters.aap.directory_conventions.analysisid_suffix;
                        thisBranchSuffix = branchStages.module(iBInd(m)).extraparameters.aap.directory_conventions.analysisid_suffix;
                        if ~strcmp(outBranchSuffix, thisBranchSuffix)
                            branchStages.module(iBInd(m)).extraparameters.aap.directory_conventions.analysisid_suffix = strcat(outBranchSuffix, thisBranchSuffix);
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
                
                branchRepeatedStages = [branchRepeatedStages branchStages.module];
                
            end % End for oB = 1 : numOutputBranches
            
            outstages.module = [outstages.module branchRepeatedStages];
            
        else
            outstages.module=extrastages.module;
        end;
        
    end; % End for loop stagenum 
    
end;
end % End Function processbranch


function outstages=checkhasrequiredfields(outstages)

% Check all required fields are present
reqflds={'index','name','epiprefix','tobecompletedfirst','extraparameters','aliasfor','remotestream','branchID'};

for stagenum=1:length(outstages.module)
    for reqfldsind=1:length(reqflds)
        if (~isfield(outstages.module(stagenum),reqflds{reqfldsind}))
            outstages.module(stagenum).(reqflds{reqfldsind})=[];
        end;
    end;
end

end
