dbstop if error

if ~exist('obj', 'var')
    global taskqueue
    obj = taskqueue;
end

if isempty(obj)
    error('No task queue found')
end

L = length([obj.pool.Jobs.ID]);

errors = false(1,L);
for ii = 1:L
    if ~isempty(obj.pool.Jobs(ii).Tasks.Error)
        errors(ii) = true;
    end
end

ei = find(errors);

%%

if ei
    for eii = ei
        aa_doprocessing_onetask(obj.pool.Jobs(eii).Tasks.InputArguments{:})
    end
else
    disp('No error found. If you are sure there was an error, then try running in localsingle instead of qsub. There may be some problems getting the job to start remotely.')
end

