% aas_makedir(aap,dirname)
% checks a directory is made - if not makes it
% Rhodri Cusack MRC CBU Cambridge Nov 2006

function [aap,resp]=aas_makedir(aap,dirname)
if (~length(dir(dirname)))
    try
        mkdir(dirname);
    catch
        aas_log(aap,1,sprintf('Problem making directory%s',dirname));
    end;
end;
