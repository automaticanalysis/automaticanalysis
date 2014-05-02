function aap = aas_configforSPM12(aap)
% T1 template
aap.directory_conventions.T1template = 'toolbox/OldNorm/T1.nii';

% TPM
modules = fieldnames(aap.tasksettings);
for i = 1:numel(modules)
    if isfield(aap.tasksettings.(modules{i}),'tpm')
        aap.tasksettings.(modules{i}).tpm = fullfile(spm('dir'), 'tpm/TPM.nii');
    end
end