function aa_test_getreport(testdir,reportdir)
    aas_makedir([],reportdir);
    projectsdirs = cellstr(spm_select('FPList',testdir,'dir','.*'));
    for d = projectsdirs(cellfun(@(p) ~isempty(spm_select('List',p,'^aap_parameters_reported.mat$')), projectsdirs))'
        aa_report_export(d{1},fullfile(reportdir,spm_file(d{1},'basename')));
    end
end