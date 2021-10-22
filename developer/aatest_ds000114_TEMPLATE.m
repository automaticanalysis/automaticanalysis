function aatest_ds000114_TEMPLATE(parameterfile, deleteprevious, wheretoprocess)

% This is an example script that can be included in the collection
% that is run using aa_test to test code submitted in a pull request.
% (This script itself isn't part of the test collection -- it's 
% only provided for documentation purposes.)
%
% There are 3 main issues when writing and using a testscript:
%
%   1) testscript/tasklist file naming and organization
%   2) location of the testscript/tasklist in the repository
%   3) parameter file customization
%
% We describe these in turn:
%
% 1) A testscript mfile must follow the naming convention:
%
%    aatest_<dataname>_<testname>.m
%   
% where "dataname" is a toplevel BIDS directory name (usually an
% OpenNeuro accession number such as "ds000114") and "testname" is
% a descriptor of the functionality being tested (examples can
% be found in $AAHOME/developer/testscripts). aa_test uses these
% strings to locate the data and results directories for the script 
% (described below).
%
% The testscript must be defined as a function taking three parameters:
% "parameterfile" is a fullpath to the aa parameterfile, "deleteprevious" 
% controls whether any existing (previous) analysis results are deleted 
% prior to running any new analysis, and "wheretoprocess" is the standard 
% processing option of localsingle, parpool, qsub, etc. The values passed 
% are determined by the options specified when aa_test is run (aa_test 
% will provide default values when necessary). See the comment block at 
% the top of aa_test.m for details. Processing of these arguments is 
% handled by the helper function aa_test_inittest. A call to this function 
% should appear as the first line in your script (see below).
%
% The function name should be the same as the filename, per the usual 
% Matlab convention.
%
% A tasklist associated with the testscript must also be provided
% with and follow the naming convention:
%
%    aatest_<dataname>_<testname>.xml
%
% i.e., the same name as the testscript except with an xml extension. 
% This is just an ordinary aa tasklist that exercises the functionality
% of interest. No special entries or modules are required.
%
% 2) The testscript and its associated tasklist must live in:
%  
%   $AAHOME/developer/testscripts
%
% otherwise aa_test won't be able to find them (the example script you are 
% currently reading lives in $AAHOME/developer -- the wrong place! -- as 
% it's only intended for documentation and not actual testing). 
%
% Note any scripts you add to your local repo for testing can be included
% in the official aa repository by opening a pull request.
%
% 3) Some customization of the aa parameter file is required for PR testing. 
% First, the following entries *must* be set in the parameter file:
%
%   aap.directory_conventions.rawdatadir
%   aap.acq_details.root
%
% Many aa users do not set these in their parameter file and instead 
% set them in the analysis script. When using aa_test, they must be
% set in the parameter file. Also, these fields are interpreted slightly
% differently by aa_test than in a standard aa analysis.
%
% Recall, aa writes analysis results to
%
%   aap.acq_details.root/aap.directory_conventions.analysisid
%
% However, when running under aa_test, the analysisid field is derived 
% from the testscript name. For example, when running the script 
% aatest_ds000456_dartel.m, results will be written to
%
%       /home/results/ds000456_dartel
%
% assuming aap.acq_details.root = '/home/results'
%
% Data used for a given script is located by appending the "dataname"  
% component of the script filename to aap.directory_conventions.rawdatadir.
% For example, when running the script aatest_ds000456_dartel, aa_test
% will look for the data in:
%
%       /home/data/PR_TESTING/ds000456
%
% assuming aap.directory_conventions.rawdatadir = '/home/data/PR_TESTING'.
% As such, rawdatadir is interpreted as the *parent* directory of the data
% directory when using aa_test (and not the data directory itself as it 
% ususally is). This "reinterpretation" provides a simple way for users
% to organize data for PR testing (a collection of BIDS directories living
% under rawdatadir).
%
% The remainder of the entries in the parameter file work as usual.
%
% Because of these customizations, you may find it convenient to create a 
% separate parameter file used only for PR testing. A fullpath to this 
% file can then be passed to aa_test. Alternatively, you can maintain
% one parameter file and just be sure to (re)set rawdatadir, root, and 
% analysisid in the script when doing a "standard" aa analysis.
%


% -------------------------------------------------------------------------
% ********************** CODE BEGINS HERE *********************************
% -------------------------------------------------------------------------

% the data and results directories are setup up by calling the convenience
% function aa_test_inittest, which also processes the options passed by
% aa_test:

% % % % % aap = aa_test_inittest(mfilename('fullpath'),deleteprevious);
aap = aa_test_inittest(mfilename('fullpath'), parameterfile, deleteprevious, wheretoprocess);

% from here down, everything works like a normal aa script. You set options
% then call aas_processBIDS with whatever parameters are appropriate
% to test the target functionality (this example just runs the motor task
% from ds000114 on one subject. Your testing should probably be more
% rigorous). 

% A report is generated if aap.directory_conventions.reportname is defined.

% -------------------------------------------------------------------------
% analysis options
% -------------------------------------------------------------------------

% % % % aap.options.wheretoprocess = wheretoprocess;

aap.options.NIFTI4D = 1;
aap.acq_details.numdummies = 4;	
aap.acq_details.input.correctEVfordummies = 1;

aap.acq_details.input.combinemultiple = 1;
aap.options.autoidentifystructural_choosefirst = 1;

% -------------------------------------------------------------------------
% BIDS input
% -------------------------------------------------------------------------

aap = aas_processBIDS(aap, [], {'finger_foot_lips'}, {'sub-01'}); 

% -------------------------------------------------------------------------
% run and report
% -------------------------------------------------------------------------

aa_doprocessing(aap);

if isfield(aap.directory_conventions,'reportname') && ~isempty(aap.directory_conventions.reportname)
    aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
end


aa_close(aap);
