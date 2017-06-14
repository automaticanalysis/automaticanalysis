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
fprintf('%d errors found\n', L)

if ei
    for eii = 1:ei
        modnum = obj.pool.Jobs(eii).Tasks.InputArguments{3};
        subjnum = obj.pool.Jobs(eii).Tasks.InputArguments{4}(1);
        prompt = sprintf('Error found in subject: %s | Module: %s\nDo you want to run this module locally y/n?\n', obj.aap.acq_details.subjects(subjnum).subjname, obj.aap.tasklist.main.module(modnum).name);
        s = input(prompt, 's');
        switch s
            case 'y'
                aa_doprocessing_onetask(obj.pool.Jobs(eii).Tasks.InputArguments{:})
            case 'n'
                disp('Skipped...')
            case 'q'
                return
        end
    end
else
    disp('No error found. If you are sure there was an error, then try running in localsingle mode instead of qsub. There may be some problems getting the job to start remotely.')
end

