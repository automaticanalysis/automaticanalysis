function aas_dumpifdeployed(aap,resp)

global aaworker;

if (isdeployed)
    try
        aapth=aaworker.parmpath;
    catch
        aapth=aaworker_getparmpath(aap,0);
    end;
    tmpvarsfn=fullfile(aapth,'tmpvars_out.mat');
    save(tmpvarsfn,'aap','resp'); 
end;