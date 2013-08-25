function aap=aas_checkspmrunning(aap)
try
    if (isempty(spm_figure('FindWin')))
        spm('fmri');
    end;
catch
    spm('fmri');
end;
end