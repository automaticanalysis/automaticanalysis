% AA module - collects garbage files after the end of an experiment
% Optional, if you don't want garbage collection to occur at each stage of
% your analysis...

function [aap,resp]=aamod_garbagecollection(aap,task)

resp='';

switch task
    case 'doit'
        
        aas_garbagecollection(aap,true);
        
    case 'checkrequirements'
        aas_log(aap,0,'Garbage needs to be taken out\n' );
end
