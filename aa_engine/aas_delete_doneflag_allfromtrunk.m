function [deps]=aas_delete_doneflag_allfromtrunk(aap,modulenum,domain)

if strcmp(domain,'study')
    deps={{'study',[]}};
else
    pth=aas_getstudypath(aap,modulenum);
    
    domaintree=aas_dependencytree_finddomain(domain,aap.directory_conventions.parallel_dependencies.study,{});
    deps=allfromtrunk(aap,pth,modulenum,domaintree,[],{});
end;

end

function [deps]=allfromtrunk(aap,pth,modulenum,domaintree,indices,deps)

if exist(pth,'dir')
    N=aas_getN_bydomain(aap,domaintree{1},indices);
    if ~isnan(N)
        if length(domaintree)>1
            for Nind=1:N
                directory=aas_getdirectory_bydomain(aap,domaintree{1},[indices Nind]);
                deps=allfromtrunk(aap,fullfile(pth,directory),modulenum,domaintree(2:end),[indices Nind],deps);
            end;
        else
            nme=aas_doneflag_getname(aap,modulenum);
            fprintf('Would remove %s\n',fullfile(pth,nme));
        end;
    end;
else
    fprintf('No directory %s so quitting tree traverse\n',pth);
end;
end