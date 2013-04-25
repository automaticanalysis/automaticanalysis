function [deps]=aas_dependencytree_allfromtrunk(aap,domain)

if strcmp(domain,'study')
    deps={{'study',[]}};
else
    domaintree=aas_dependencytree_finddomain(domain,aap.directory_conventions.parallel_dependencies.study,{});
    deps=allfromtrunk(aap,domaintree,[],{});
end;

end

function [deps]=allfromtrunk(aap,domaintree,indices,deps)
N=aas_getN_bydomain(aap,domaintree{1},indices);
if ~isnan(N)
    if length(domaintree)>1
        for Nind=1:N
                deps=allfromtrunk(aap,domaintree(2:end),[indices Nind],deps);
        end;
    else
        newdeps=cell(N,1);
        for Nind=1:N
            newdeps{Nind}={domaintree{1} [indices Nind]};
        end;
        deps=[deps;newdeps];
    end;
end;
end