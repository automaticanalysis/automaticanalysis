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
%   fashion. Other options include parpool (PCT) and qsub (cluster)
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
argParse.addParameter('deleteprevious', false, @(x) islogical(x) || isnumeric(x));
argParse.addParameter('haltonerror', false, @(x) islogical(x) || isnumeric(x));
argParse.addParameter('wheretoprocess','localsingle', @ischar);
argParse.addParameter('parameterfile','', @ischar);
argParse.parse(varargin{:});

% parse glob negation tilde

globflag = 0;
glob = argParse.Results.glob;
if ~isempty(glob)
    globflag = 1;
    if startsWith(glob,'~')
        globflag = -1;
        glob = glob(2:end);
    end
end

% logging

if exist(fullfile(pwd,'aa_test.log'),'file') 
    ow = input('aa_test.log exists. Overwrite?(Y/[N])> ','s');
    if (isempty(ow) || ow == 'N' || ow == 'n'); return; end
end

fid = fopen(fullfile(pwd,'aa_test.log'),'w');

if (fid < 0)
    error('Cannot open aa_test.log for writing');
end
   
% get a list of testscripts
aa = aaClass('nopath','nogreet');
aadir = aa.Path;
testdir = fullfile(aadir,'developer','testscripts');

tests = dir(fullfile(testdir,'*.m'));
testnames = {tests.name};

savedir = pwd;
cd(testdir);

for fname = testnames
    if globflag<0 && contains(fname{1},glob); continue; end
    if globflag>0 && ~contains(fname{1},glob); continue; end
    fprintf('\n\nRunning %s...\n', fname{1});
    runit(fname{1}, argParse.Results.parameterfile, argParse.Results.deleteprevious, argParse.Results.haltonerror, argParse.Results.wheretoprocess, fid);
end

fprintf('------------------\ntests finished.\n')
cd(savedir);
fclose(fid);

end

% -------------------------------------------------------------------------
function runit(scriptname, parameterfile, deleteprevious, haltonerror, wheretoprocess, fid)
% -------------------------------------------------------------------------

func = str2func(strrep(scriptname,'.m',''));

if  haltonerror
    
    % run script normally;
    % function halts on aa error
    
    func(parameterfile, deleteprevious, wheretoprocess);
    
    % if encapsulated returns, this script passed
    
    fprintf('\npass - %s\n',scriptname);
    fprintf(fid,'\npass - %s\n',scriptname);

else
    
    % run script inside a try-catch;
    % function flags any error and returns
    % optionally dropping into qsub_debug
    
    try
        
        func(parameterfile, deleteprevious, wheretoprocess);
        fprintf('\npass - %s\n',scriptname);
        fprintf(fid,'\npass - %s\n',scriptname);
        
    catch err
        
        fprintf('FAIL - %s\n',scriptname);
        fprintf(fid,'FAIL - %s\n',scriptname);
        
        if ~strcmp(wheretoprocess,'localsingle')
            
            try
                aaq_qsub_debug;
            catch
                fprintf('aaq_qsub_debug: %s. See err variable.\n',err.message);
                fprintf(fid,'aaq_qsub_debug: %s. See err variable.\n',err.message);
            end
            
            % catch cases when aaq_qsub_debug quietly 
            % returns (e.g., failed to submitjob)
            
            keyboard;
            
        end
        
    end

end

end