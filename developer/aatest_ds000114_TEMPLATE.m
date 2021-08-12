function aatest_ds000114_TEMPLATE(deleteprevious,wheretoprocess)
% developer test script
%
% this is a minimal script used to document how to setup
% a script for PR testing using the new version of aa_test
%
% data used:  ds000114 (OpenNeuro)
%
% (the accompanying tasklist aatest_ds000114_TEMPLATE.xml
% simply converts the structural and exits. Actual PR testing
% scripts will of course include more extensive functionality)
%

% -------------------------------------------------------------------------
% init
% -------------------------------------------------------------------------

aap = aa_test_inittest(mfilename('fullpath'),deleteprevious);

% Scripts for PR testing load parameterset ~/.aa/aap_parameters_user.xml,
% which may include a site-specific or the default parameterset. You should
% customize this or the site-specific parameterset for how you want to do
% PR testing.
%
% However, two settings comes from the script filename, and we assume a 
% testscript filename has format: 
% aatest_<BIDS directory name>_<test name>.m:
%
% 1) the BIDS directory name 
% 2) the test name describing the nature of the test
%
% We assume the BIDS data for PR testing lives in 'aa_demo' folder listed
% in aap.directory_conventions.rawdatadir
% (i.e. /<some path>/aa_demo/<BIDS directory name>).
%
% The result folder as <BIDS directory name>_<test name> will be created
% under aap.acq_details.root
% 
% For example running the scripts
%
%   aatest_ds000123_quicktest.m
%   aatest_ds000123_fulltest.m
%   aatest_ds000456_dartel.m
%
% uses datasets in
%
%   /home/data/aa_demo/ds000123
%   /home/data/aa_demo/ds000123
%   /home/data/aa_demo/ds000456
%
% assuming aap.directory_conventions.rawdatadir contains '/home/data/aa_demo'
%
% additionally, the result will be placed in
%
%   /home/results/ds000123_quicktest
%   /home/results/ds000123_fulltest
%   /home/results/ds000456_dartel
%
% assuming aap.acq_details.root = '/home/results'
%
% note if ~/.aa/aap_parameters_user.xml doesn't exist, aa will try to
% create one for you when the script runs, but it almost certainly
% require further customization to work. The take-home is just create 
% a proper parameter file in ~/.aa/aap_parameters_user.xml before 
% starting your PR testing.
%
% the tasklist .xml files are assumed to live in the same directory with the same
% basename as the scripte
%
% from here down, everything should work normally...

% -------------------------------------------------------------------------
% analysis options
% -------------------------------------------------------------------------

aap.options.wheretoprocess = wheretoprocess;

aap.options.NIFTI4D = 1;
aap.acq_details.numdummies = 4;	
aap.acq_details.input.correctEVfordummies = 1;

aap.acq_details.input.combinemultiple = 1;
aap.options.autoidentifystructural_choosefirst = 1;

% -------------------------------------------------------------------------
% BIDS
% -------------------------------------------------------------------------

% the script developer needs to decide what task, sessions, and subjects 
% are required to test the target functionality 

aap = aas_processBIDS(aap, [], {'finger_foot_lips'}, {'sub-01'}); 

% -------------------------------------------------------------------------
% run
% -------------------------------------------------------------------------

aa_doprocessing(aap);

% if directory_conventions.reportname is undefined, skip reporting

if isfield(aap.directory_conventions,'reportname') && ~isempty(aap.directory_conventions.reportname)
    aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
end


aa_close(aap);
