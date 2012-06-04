function [domaintree]=aas_dependencytree_finddomain(domain,tree,path)

fn=fieldnames(tree);

domaintree={};
for fnind=1:length(fn)
    if strcmp(fn{fnind},domain)
        domaintree={path{:} domain};
        return;
    else
        domaintree=aas_dependencytree_finddomain(domain,tree.(fn{fnind}),{path{:} fn{fnind}});
        if ~isempty(domaintree)
            return;
        end;
    end;
end;
end