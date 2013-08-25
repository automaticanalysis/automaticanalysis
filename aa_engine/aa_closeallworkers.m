% aa parallel
% Kills all workers. They'll decay anyway, but sometimes nice to clean them
% out
%

function aa_closeallworkers()
global aaparallel

clear aap;


if (~isempty(aaparallel))
    aaprocesskey=aaparallel.processkey;
    aap.internal.parallel.processkey=aaprocesskey;
    pth=aaworker_getparmpath(aap,0,true);
    [subpth nme ext]=fileparts(pth);
    fn=dir(fullfile(subpth,['aaworker' num2str(aaprocesskey) '*']));
    for i=1:length(fn)
        if (fn(i).name(1)~='.')
            fprintf('Searching for processes associated with %s\n',fn(i).name);
            workerprocesskey=str2num(fn(i).name(9:end));
            aas_sackworker(aaparallel,workerprocesskey);
        end;
    end;
    aaparallel.workerlist=[];
end;