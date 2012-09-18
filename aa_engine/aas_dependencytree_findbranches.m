% Traverse 'tree' (a structure), picking out all instances at each level -
% e.g., every subject or session
%  If optional parameter '
% First level to be fully traversed is all those fields of tree
%  So, for example, if structure is:
%   tree.subjects.sessions
%  then all subjects and all sessions of those subjects will be returned
%  The indices should correspond to the parts of the tree not provided (i.e. no
%  indices in the above example)
%
% If you provide a modulenum then it assumes you only want to return
% branches for which the subdirectory exists
%  i.e., where there might be a doneflag

function [deps]=aas_dependencytree_findbranches(aap,tree,indices,onlyreportdomain,modulenum,deps)
if ~exist('modulenum','var')
    modulenum=[];
end;
if ~exist('deps','var')
    deps={};
end;
if ~exist('onlyreportdomain','var')
    onlyreportdomain=[];
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
            N=aas_getN_bydomain(aap,fn{fnind},indices);
            range=1:N;
        end;
        for branchind=range
            traverse=true;
            if exist('modulenum','var') && ~isempty(modulenum)
                pth=aas_getpath_bydomain(aap,fn{fnind},[indices branchind],aap.directory_conventions.remotefilesystem,modulenum);
                if ~exist(pth,'dir')
%                    fprintf('Quitting traverse as directory does not exist:\n %s\n',pth);
                    traverse=false;
                end;
            end;
            if traverse
                
                % Add this
                if isempty(onlyreportdomain) || strcmp(onlyreportdomain,fn{fnind})
                    deps{end+1}={fn{fnind} [indices branchind]};
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
                    deps=aas_dependencytree_findbranches(aap,treecell,[indices branchind],onlyreportdomain,modulenum,deps);
                end;
            end;
        end;
    end;
end;

end