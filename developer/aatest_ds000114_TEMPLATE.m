
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

clear;
aa_ver5;

% scripts for PR testing load parameter file ~/.aa/aap_parameters_user.xml
% you should customize this for how you want to do PR testing
%
% However, two settings comes from the script filename:
%
% 1) the BIDS directory name (appended to aap.directory_conventions.rawdatadir) 
% 2) the results directory name (created under aap.acq_details.root --
% i.e., aap.directory_conventions.analysisid = script_name (w/o extension)
% 
% note assumption 1 means you must keep all of your PR testing data udner
% aap.directory_conventions.rawdatadir specified in your "PR testing"
% parameters script (~/.aa/aap_parameters_user.xml). 
%
% For example running the scripts
%
%   aatest_ds000123_quicktest
%   aatest_ds000123_fulltest
%   aatest_ds000456_dartel
%
%   /home/data/PR_TESTING/ds000123
%   /home/data/PR_TESTING/ds000123
%   /home/data/PR_TESTING/ds000456
%
% assuming aap.directory_conventions.rawdatadir = /home/data/PR_TESTING
% in ~/.aa/aap_parameters_user.xml
%
% additionally, the result will be placed in
%
%   /home/results/PR_RESULTS/ds000123_quicktest
%   /home/results/PR_RESULTS/ds000123_fulltest
%   /home/results/PR_RESULTS/ds000456_dartel
%
% assuming aap.acq_details.root = '/home/results/PR_RESULTS';
% in ~/.aa/aap_parameters_user.xml
%

% note if ~/.aa/aap_parameters_user.xml doesn't exist, aa will try to
% create one for you when the script runs, but it almost certainly
% require further customization to work. The take-home is just create 
% a proper parameter file in ~/.aa/aap_parameters_user.xml before 
% starting your PR testing.

% the .xml is assumed to live in the same directory as this script

aap = aarecipe([mfilename('fullpath') '.xml']);

% -------------------------------------------------------------------------
% directory setup
% -------------------------------------------------------------------------

% we assume a testscript filename has format: aatest_foo_bar.m
% we assume data lives in aap.directory_conventions.rawdatadir/foo
% we assume results should go in aap.acq_details.root/foo_bar

% aap.acq_details.root =  <<<< set this in ~/.aa/aap_parameters_user.xml! >

temp = split(mfilename,'_');
aap.directory_conventions.analysisid = [ temp{2} '_' temp{3} ];

fprintf('Saving results in: %s/%s\n', aap.acq_details.root, aap.directory_conventions.analysisid);

% NB: meaning of aap.directory_conventions.rawdatair
%
% usually its the directory where the data lives -- here its one level up
% because we get the actual directory from the script name. 
%
% For example, usually we would have
%
% aap.directory_conventions.rawdatadir = '/path/to/ds000123'
% (where ds000123 is a BIDS directory)
%
% but for PR testing, do: 
%
%   aap.directory_conventions.rawdatadir = '/path/to'
%
% because the script needs to append ds000123 here
%
% this convention is the only way we can use a common set of
% test scripts without requiring everybody to put data in
% some (inconvenient) hard coded named directory

aap.directory_conventions.rawdatadir = fullfile(aap.directory_conventions.rawdatadir,temp{2});

% be sure to set this up in ~/.aa/aap_parameters_user.xml
% aap.options.wheretoprocess = ... 

% from here down, everything should work normally...

% -------------------------------------------------------------------------
% analysis options
% -------------------------------------------------------------------------

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
