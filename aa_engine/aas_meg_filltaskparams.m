function [aap]=aas_meg_filltaskparams(aap, task, I)

fields=fieldnames(I);
for f=1:length(fields)
    if ~isfield(aap.MEG.TaskSettings.(task),fields{f})
        % use default
        aap.MEG.TaskSettings.(task).(fields{f})=I.(fields{f});
    end
end
