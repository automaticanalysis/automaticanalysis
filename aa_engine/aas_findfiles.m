%function [aap,resp]=aas_findfiles(aap,searchpath,minfiles,maxfiles)
% if minfiles is not specified it defaults to 1
% if maxfiles is not specified it defaults to minfiles
% checks a directory is made - if not makes it
% Rhodri Cusack MRC CBU Cambridge Nov 2006

function [aap,resp]=aas_findfiles(aap,searchpath,minfiles,maxfiles)
cmd=['ls ' searchpath];
[s fles]=aas_shell(cmd);
if (nargin<3)
    minfiles=1;
end;
if (nargin<4) 
    maxfiles=minfiles;
end;
nfles=size(fles,1);
if (nfles<minfiles)
    aas_log(aap,1,sprintf('Only %d files were found like \n%s\nwhich is not enough because the minimum is %d',nfles,searchpath,minfiles));
end;
if (nfles>maxfiles)
    aas_log(aap,1,sprintf('Found %d files like \n%s\nwhich is too many because the maxmimum is %d',nfles,searchpath,maxfiles));
end;
[pth nme ext]=fileparts(searchpath);
resp=fles;