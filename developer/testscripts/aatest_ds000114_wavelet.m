function aatest_ds000114_wavelet(parameterfile, deleteprevious, wheretoprocess)

% This script tests aamod_waveletdespike on a single subject using 
% the ds000114 motor task. The data should be downloaded prior
% to running the script

% see aatest_ds00114_TEMPLATE.m in $AAHOME/developer for help on
% writing and using a test script

% -------------------------------------------------------------------------
% init
% -------------------------------------------------------------------------

aap = aa_test_inittest(mfilename('fullpath'), parameterfile, deleteprevious, wheretoprocess);

% -------------------------------------------------------------------------
% analysis options
% -------------------------------------------------------------------------

aap.options.autoidentifystructural_choosefirst = 1;
aap.options.autoidentifystructural_chooselast = 0;

aap.options.NIFTI4D = 1;

aap.acq_details.numdummies = 4;	
aap.acq_details.input.correctEVfordummies = 1;

aap.tasksettings.aamod_segment8.writenormimg = 0;
aap.tasksettings.aamod_segment8.samp = 2;
aap.tasksettings.aamod_smooth.FWHM = 5;

% -------------------------------------------------------------------------
% BIDS
% -------------------------------------------------------------------------

aap.acq_details.input.combinemultiple = true;
aap = aas_processBIDS(aap,[],{'finger_foot_lips'},{'sub-01'});

% -------------------------------------------------------------------------
% modeling - contrast specification
% -------------------------------------------------------------------------

aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts_00001', '*', 'sameforallsessions', [1 0 0], 'Finger','T');
aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts_00001', '*', 'sameforallsessions', [0 1 0], 'Foot', 'T');
aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts_00001', '*', 'sameforallsessions', [0 0 1], 'Lips', 'T');

% -------------------------------------------------------------------------
% run
% -------------------------------------------------------------------------

aa_doprocessing(aap);

% if directory_conventions.reportname is undefined, skip reporting

if isfield(aap.directory_conventions,'reportname') && ~isempty(aap.directory_conventions.reportname)
    aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
end


aa_close(aap);


