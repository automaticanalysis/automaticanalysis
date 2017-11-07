% function [deps desc]=aas_dependencytree_findbranches(aap,tree,indices,onlyreportdomain,modulenum,deps,requiredancestor,hasancestory,onlyrequiredancestor)
%
% Traverse 'tree' (a structure), picking out all instances at each level -
% e.g., every subject or session
% 
% First level to be fully traversed is all those fields of tree
%  So, for example, if structure is:
%   tree.subjects.sessions
%  then all subjects and all sessions of those subjects will be returned
%  The indices should correspond to the parts of the tree not provided 
%
% onlyreportdomain can be a string, specifying the single domain to be reported 
%
% If you provide a modulenum then it assumes you only want to return
% branches for which the subdirectory exists
%  i.e., where there might be a doneflag already written
%
% If two outputs, the 'desc' is filled with a text description of the
% dependencies. This is useful for debugging
%
% Set hasancestory=false and requiredancestor='adomain' to return only
% dependencies at or beneath this ancestor.
%
% e.g., [deps desc]=aas_dependencytree_findbranches(aap,aap.directory_conventions.parallel_dependencies.study.subject,[1])
%  desc='session[ 1 1 ] session[ 1 2 ] session[ 1 3 ] session[ 1 4 ]'
% e.g., [deps desc]=aas_dependencytree_findbranches(aap,aap.directory_conventions.parallel_dependencies.study,[])
%  desc='subject[ 1 ] session[ 1 1 ] session[ 1 2 ] session[ 1 3 ] session[ 1 4 ] subject[ 2 ] session[ 2 1 ] session[ 2 2 ] session[ 2 3 ] session[ 2 4 ]';
% e.g., [deps desc]=aas_dependencytree_findbranches(aap,aap.directory_conventions.parallel_dependencies.study,[],'subject')
%  ...which had onlyreportdomain='subject', so get...
%  desc='subject[ 1 ] subject[ 2 ] subject[ 3 ] subject[ 4 ] subject[ 5 ]';


function [deps desc]=aas_dependencytree_findbranches(aap,tree,indices,onlyreportdomain,modulenum,deps,requiredancestor,hasancestory)
if ~exist('modulenum','var')
    modulenum=[];
end;
if ~exist('deps','var')
    deps={};
end;
if ~exist('requiredancestor','var')
    requiredancestor=[];
end;
if ~exist('hasancestory','var')
    hasancestory=true;
end;
if ~exist('onlyreportdomain','var')
    onlyreportdomain=[];
end;
if ~exist('onlyrequiredancestor','var')
    onlyrequiredancestor=false;
end;

% Can provide tree as cell array with 3 elements, first of which is a tree,
% second a list of filter domains, and third indices for those domains.
if iscell(tree)
    filter_domainlist=tree{2};
    filter_indices=tree{3};
    tree=tree{1};
else
    filter_domainlist=[];
    filter_indices=[];
end;


% First we need to go through every node
if isstruct(tree)
    fn=fieldnames(tree);
    for fnind=1:length(fn)
        % If we've been told to pick a particular branch for a domain at
        % this level, do it.
        if ~isempty(filter_domainlist) && strcmp(fn{fnind},filter_domainlist{1})
            range=filter_indices(1);
            filter_domainlist=filter_domainlist(2:end);
            filter_indices=filter_indices(2:end);
        else
            % otherwise, provide every branch
            [N, range]=aas_getN_bydomain(aap,fn{fnind},indices);
        end;
        for branchind=range
            traverse=true;
            if exist('modulenum','var') && ~isempty(modulenum)
                pth=aas_getpath_bydomain(aap,fn{fnind},[indices branchind],aap.directory_conventions.remotefilesystem,modulenum);
                if ~exist(pth,'dir')
                    % aas_log(aap,false,sprintf('INFO: Quitting traverse as directory does not exist:\n\t%s',pth));
                    traverse=false;
                end
            end
            if traverse
                
                % Have we reached the required ancestor? 
                if strcmp(fn{fnind},requiredancestor)
                    hasancestory=true;
                end;
                
                % Add this
                if hasancestory
                    if isempty(onlyreportdomain) || strcmp(onlyreportdomain,fn{fnind}) 
                        deps{end+1}={fn{fnind} [indices branchind]};
                    end;
                end;
                
                % And go through its branches, making sure not to traverse
                % lower than a requested domain if specified in
                % onlyreportdomain
                if isempty(onlyreportdomain) || ~strcmp(onlyreportdomain,fn{fnind})
                    if ~isempty(filter_domainlist)
                        treecell=tree.(fn{fnind});
                    else
                        treecell={tree.(fn{fnind}) filter_domainlist filter_indices};
                    end;
                    
                    
                        deps=aas_dependencytree_findbranches(aap,treecell,[indices branchind],onlyreportdomain,modulenum,deps,requiredancestor,hasancestory);
                    
                end;
            end;
        end;
    end;
end;

%% And report
if nargout==2
    desc='';
    for depind=1:length(deps)
        desc=[desc deps{depind}{1} '[ ' sprintf('%d ',deps{depind}{2}) '] '];
    end;
end;
end