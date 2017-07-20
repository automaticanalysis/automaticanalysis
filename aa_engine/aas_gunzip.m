% Gunzip a .gz file
% function [aap fn]=aas_gunzip(aap,fn)

function [aap fn]=aas_gunzip(aap,fn)

[pth nme ext]=fileparts(fn);
if strcmp(ext,'.gz')
    if exist(fn,'file')
        [s w]=aas_shell(['gunzip -f ' fn]);
        if s
            aas_log(aap,true,sprintf('Error unpacking %s',fn));
        end;
        fn=fullfile(pth,nme);
    else
        aas_log(aap,false, sprintf('aas_gunzip did''t find %s, already unpacked?',fn));
    end;
end;