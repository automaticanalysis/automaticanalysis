function [deps]=aas_dependencytree_findbranches(aap,tree,deps,indices)
% First we need to go through every node
if isstruct(tree)
    fn=fieldnames(tree);
    for fnind=1:length(fn)
        % And follow every branch from that node
        N=aas_getN_bydomain(aap,fn{fnind},indices);
        for branchind=1:N
            % Add this
            deps{end+1}={fn{fnind} [indices branchind]};
            % And go through its branches
            deps=aas_dependencytree_findbranches(aap,tree.(fn{fnind}),deps,[indices branchind]);
        end;
    end
else
    N=nan;
end;
end