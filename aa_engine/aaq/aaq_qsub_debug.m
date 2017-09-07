function aaq_qsub_debug()

global taskqueue
obj = taskqueue;

dbstop if error

if isempty(obj)
    dbclear all
    error('No task queue found. This function can only be used after receiving an error running AA in "qsub" mode.')
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
    fprintf('%d errors found\n', length(ei))
    for eii = ei
        modnum = obj.pool.Jobs(eii).Tasks.InputArguments{3};
        if obj.pool.Jobs(eii).Tasks.InputArguments{4}
            subjnum = obj.pool.Jobs(eii).Tasks.InputArguments{4}(1);
            prompt = sprintf('Error found in subject: %s | Module: %s\nDo you want to run this module locally y/n?\n',...
                obj.aap.acq_details.subjects(subjnum).subjname, obj.aap.tasklist.main.module(modnum).name);
        else
            prompt = sprintf('Error found in Session Level Module: %s\nDo you want to run this module locally y/n?\n',...
                obj.aap.tasklist.main.module(modnum).name);
        end
        
        s = input(prompt, 's');
        switch s
            case 'y'
                aa_doprocessing_onetask(obj.pool.Jobs(eii).Tasks.InputArguments{:})
            case 'n'
                disp('Skipped...')
            case 'q'
                dbclear all
                disp('Clearing breakpoints and exiting debug mode!')
        end
    end
else
    disp('No error found. If you are sure there was an error, then try running in localsingle mode instead of qsub. There may be some problems getting the job to start remotely.')
end

dbclear all