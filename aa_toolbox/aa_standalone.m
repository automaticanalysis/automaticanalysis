function aap = aa_standalone(fname_aa, fname_tasklist)
aap=aarecipe('aap_parameters_defaults_CBSU.xml',fname_tasklist);
run(fname_aa)

%% DO ANALYSIS
aa_doprocessing(aap);
% aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
