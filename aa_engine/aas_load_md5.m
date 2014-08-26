function [aap md5 datecheck isarchived]=aas_load_md5(aap,fid,streamname)
isarchived=false;
if (ischar(fid))
    lne=fid;
else
    lne=fgetl(fid);
    
    % Has this file been archived?
    if strcmp(lne(1:11),'ARCHIVED TO')
        lne=fgetl(fid);
        isarchived=true;
    end;
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

