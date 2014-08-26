function aas_checkspmrunning(aap,force)
if strcmp(aap.options.wheretoprocess,'localsingle') || ((nargin > 1) && force)
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