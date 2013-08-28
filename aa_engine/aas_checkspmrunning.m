function aap=aas_checkspmrunning(aap)
try
    if (isempty(spm_figure('FindWin')))
        spm('fmri');
    end;
catch
    spm('fmri');
end;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'extrafunctions','spm_mods'),'-begin');
end