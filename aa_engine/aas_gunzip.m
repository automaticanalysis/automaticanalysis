function [aap fn]=aas_gunzip(aap,fn)

[pth nme ext]=fileparts(fn);
if strcmp(ext,'.gz')
    [s w]=aas_shell(['gunzip -f ' fn]);
    if s
        aas_log(aap,true,sprintf('Error unpacking %s',fn));
    end;
    fn=fullfile(pth,nme);
end;
