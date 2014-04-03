function [aap md5 datecheck]=loadmd5(aap,fid,streamname)
if (ischar(fid))
    lne=fid;
else
    lne=fgetl(fid);
end;
pos=find(lne==9);
if (length(lne)<3 || ~strcmp(lne(1:3),'MD5') || isempty(pos))
    aas_log(aap,true,sprintf('MD5 in file %s corrupted',streamname));
end;
if (length(pos)==1)
    md5=deblank(lne(pos+1:end));
    datecheck='';
else
    md5=deblank(lne(pos(1)+1:pos(2)-1));
    datecheck=deblank(lne(pos(2)+1:end));
end;