% aa parallel
% Kills all workers. They'll decay anyway, but sometimes nice to clean them
% out
%

function aa_closeallworkers
global aaparallel

clear aap;


if (~isempty(aaparallel))
    aaprocesskey=aaparallel.processkey;
    aap.internal.parallel.processkey=aaprocesskey;    
    subpth=fileparts(aaworker_getparmpath(aap,0,true));
    fn=dir(fullfile(subpth,['aaworker' num2str(aaprocesskey) '*']));
    for i=1:length(fn)
        if fn(i).name(1)~='.'
            workerid = fn(i).name(9:end);
            fprintf('Searching for processes with ID = %s\n',workerid);
            aas_sackworker(aaparallel,workerid);
        end;
    end;
    aaparallel.workerlist=[];
end;