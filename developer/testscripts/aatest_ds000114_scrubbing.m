function aatest_ds000114_scrubbing(parameterfile, deleteprevious, wheretoprocess)

% This script tests aa frame censoring pipeline using 4
% subjects in the the ds000114 motor task. The data should be 
% downloaded prior to running the script

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

aap.tasksettings.aamod_scrub_epi.scrub_criteria = 'metric_data.DVARS > metric_thresholds.DVARS.fivepercent';

aap.tasksettings.aamod_firstlevel_model.includemovementpars = 0;
aap.tasksettings.aamod_firstlevel_model.includespikes= 1;

aap.tasksettings.aamod_firstlevel_threshold.description = '5% DVARS, 0.001 UNC';

% -------------------------------------------------------------------------
% BIDS
% -------------------------------------------------------------------------

aap.acq_details.input.combinemultiple = true;
aap = aas_processBIDS(aap,[],{'finger_foot_lips'},{'sub-01','sub-02','sub-03','sub-04'});

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


