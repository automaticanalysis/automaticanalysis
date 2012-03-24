% AA apply remap table to a stream
% Rhodri Cusack MRC CBU, Cambridge Dec 2010

function newstreamname=aas_remapstreamname(aap,oldstreamname,isinput)

try
    if (isinput)
        remap=aap.tasklist.currenttask.inputstreams.remap;
    else
        remap=aap.tasklist.currenttask.outputstreams.remap;
    end;
    
catch
    newstreamname=oldstreamname;
    remap=[];
end;

for remapind=1:length(remap)
    if (strcmp(oldstreamname,deblank(remap(remapind).from)))
        newstreamname=remap(remapind.to);
        break;
    end;
end;

