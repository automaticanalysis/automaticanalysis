function aa_test(varargin)
%
% run $AAHOME/developer testscripts
%
% inputs (optional) - specified as keyword,value pairs
%
% 'parameterfile',parameterfile -
%   == (fullpath) to parameter file. If not specified, aa_test will
%   use ~/.aa/aap_parameters_user.xml. NB: the parameter file used
%   for aa_test is slightly different than the standard parameter file.
%   Briefly, aap.directory_conventions.rawdatadir and
%   aap.acq_details.root *must* be defined in the parameter file (many
%   aa users prefer to define these in the analysis script, not the
%   parameter file) and rawdatadir is interpreted as the *parent*
%   directory of the test data directory (or directories), not the data
%   directory itself. See the template test script provided in
%   $AAHOME/developer/aatest_ds000114_TEMPLATE.m for an expanded
%   discussion (including the motivation for these changes).
%
% 'glob', glob -
%   == (limited) glob search string used to restrict which
%   testscripts are run ("limited" in the sense that we only
%   recognize a string literal and a tilde-negated string literal.
%   See example usage below). If not specified, all testscripts
%   are run.
%
% 'deleteprevious', deleteprevious -
%   == true, testing will delete any previous results to force
%   re-execution of the entire analysis (default: true)
%
% 'haltonerror', haltonerror -
%   == true, testing will halt (crash) if any errors, otherwise
%   testing continues with the next script (default: false)
%
% 'wheretoprocess', wheretoprocess-
%   == 'localsingle' (default), jobs will run locally in a serial
%   fashion. Other options include 'parpool' and 'batch' (cluster)
%
% Example Usage
%
%   aa_test;                      - run all scripts, use default settings
%   aa_test('haltonerror',true)   - run all scripts, halt on error
%
%   aa_test('parameterfile','/path/to/myPRtestingparameterfile.xml');
%
% example glob usage
%
%   aa_test('glob','aatest_ds000114_fmri') - only run aatest_ds000114_fmri.m
%   aa_test('glob','~ds002737') - run all scripts EXCEPT those with
%                          "ds002737" in the name (note leading tilde)
%
% Notes
%
% test results (pass/fail) are logged to aa_test.log in working dir.
% If this file exists, you'll be asked before it is overwritten
%
% Expected usage is: 1) do 'aa_test' to run all scripts using
% defaults, 2) check the log, then 3) run aa_test('glob','foo','haltonerror'
% ,true) or aa_test('glob','foo','haltonerror',true,'wheretoprocess','qsub')
% to run the jobs on a cluster and re-run script "foo" that failed in order
% to drop into the debugger.
%
% test scripts are assumed to live in $AAHOME/developer/testscripts
% These are designed for efficient PR testing and not necessarily
% helpful for learning aa. See $AAHOME/examples for the latter.
%
% Revision History
%
% 10/2021 [MSJ] - add parameterfile parameter
% summer/2021 [MSJ] - newish (derived from aatest)
%

argParse = inputParser;
argParse.addParameter('glob','', @ischar);
argParse.addParameter('deleteprevious', true, @(x) islogical(x) || isnumeric(x));
argParse.addParameter('haltonerror', false, @(x) islogical(x) || isnumeric(x));
argParse.addParameter('wheretoprocess','localsingle', @ischar);
argParse.addParameter('parameterfile','', @ischar);
argParse.parse(varargin{:});

thisFolder = fileparts(mfilename("fullpath"));

% logging
logfile = fullfile(thisFolder,'aa_test.log');
if exist(logfile,'file')
    ow = input('aa_test.log exists. Overwrite?(Y/[N])> ','s');
    if (isempty(ow) || ow == 'N' || ow == 'n')
        return
    else
        % Delete, because do not want to append.
        delete(logfile)
    end
end

%% Get a list of all tests

% Get Use case tests
testUseCases.pass_inputargs('set', argParse.Results);
suiteUseCases = matlab.unittest.TestSuite.fromClass(?testUseCases);

testsFolder = fullfile(thisFolder, 'tests');
suiteTests = matlab.unittest.TestSuite.fromFolder(testsFolder);

suite = [suiteUseCases, suiteTests];

%% Apply glob-based filtering
glob = argParse.Results.glob;
if ~isempty(glob)
    globflag = 1;
    % parse glob negation tilde
    if startsWith(glob,'~')
        globflag = -1;
        glob = glob(2:end);
    end
    constr = matlab.unittest.constraints.ContainsSubstring(glob);
    do_select = matlab.unittest.selectors.HasName(constr);
    if globflag < 0
        suite = suite.selectIf(~do_select);
    elseif globflag > 0
        suite = suite.selectIf(do_select);
    end
end

% PLACEHOLDER: Here, apply e.g. Tag based filtering

%% Run tests
runner = matlab.unittest.TestRunner.withTextOutput;
file_plugin = matlab.unittest.plugins.ToFile(logfile);
tap_plugin = matlab.unittest.plugins.TAPPlugin.producingOriginalFormat(file_plugin);
runner.addPlugin(tap_plugin);
results = runner.run(suite);

if argParse.Results.haltonerror
    % The unittest framework catches errors during tests
    % Here, throw an error is any test failed. A.o. to notify the
    % Continuous Integration of failure.
    assertSuccess(results);
end

end
