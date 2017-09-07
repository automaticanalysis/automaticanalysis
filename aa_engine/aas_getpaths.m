
% Recursively unpack directory structure
function [aap paths dtes szes]=aas_getpaths(aap,filelist,pth)
paths={}; dtes={}; szes=[];
for i=1:length(filelist)
    if (filelist(i).isdir)
        [aap paths_app dtes_app szes_app]=aas_getpaths(aap,filelist(i).subdir,fullfile(pth,filelist(i).name));
        paths=[paths paths_app];
        dtes=[dtes dtes_app];
        szes=[szes szes_app];
    else
        paths=[paths fullfile(pth,filelist(i).name)];
        dtes=[dtes filelist(i).date];
        szes=[szes filelist(i).bytes];
    end;
end;
end