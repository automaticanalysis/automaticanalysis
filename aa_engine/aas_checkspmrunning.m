function aap=aas_checkspmrunning(aap)
if strcmp(aap.options.wheretoprocess,'localsingle')
    try
        if (isempty(spm_figure('FindWin')))
            spm('fmri');
        end;
    catch
        spm('fmri');
    end;
end
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'extrafunctions','spm_mods'),'-begin');
end