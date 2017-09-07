% function [aap dirlist]=aas_recursivedir(aap,fn)


function [aap dirlist]=aas_recursivedir(aap,fn,filefilter)

if (~exist(fn,'file'))
    dirlist=[];
else
    dirlist=dir(fn);
    
    finddirs=find([dirlist.isdir]);
    
    keep=true(length(dirlist),1);
    
    for subind=1:length(finddirs)
        nme=dirlist(finddirs(subind)).name;
        if (~strcmp(nme,'.') && ~strcmp(nme,'..'))
            [aap dirlist_sub]=aas_recursivedir(aap,fullfile(fn,nme));
            dirlist(finddirs(subind)).subdir=dirlist_sub;
        else
            keep(finddirs(subind))=false;
        end;
    end;
    
    dirlist=dirlist(keep);
end
end