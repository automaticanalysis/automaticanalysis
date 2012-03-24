function [aap waserror respstruct]=aas_runpython_getstruct(aap,pythoncode,continueif404)
if (~exist('continueif404','var'))
    continueif404=false;
end;

retrydelays=[0.5 1 2 4 8 16 32 64];
for delay=retrydelays
    try
        [aap waserror respxml]=aas_runpython(aap,pythoncode,continueif404);
        [aap respstruct]=aas_xmltostruct(aap,respxml);
        break;
    catch
        pause(delay);
    end;
end;
        
