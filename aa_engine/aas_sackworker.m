function [aap]=aas_sackworker(aap,workerprocesskey);
% Now recursive kill of child jobs as well

global aaparallel

% Get directory for worker directories
parmpth=aaworker_getparmpath(aap,workerprocesskey,true);

%
try
    aaworker=xml_read(fullfile(parmpth,'aaworker.xml'));
    try
        aaworker.pid;
    catch
    end;

    % Now kill all of the pids
    %   use kill recursive now!
    aas_recursivekill(aap,aaworker.hostname,aaworker.pid);
    % And remove directory
    try
        rmdir(parmpth,'s');
    catch
    end;

catch
end;


% And remove
try
    aaparallel.workerlist(aaparallel.workerlist==workerprocesskey)=[];
catch
end;
