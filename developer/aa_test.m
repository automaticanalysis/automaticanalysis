function aa_test(glob,haltonerror,useqsubdebug)
%
% run /developer testscripts (new version -- replaces aatest)
%
% inputs (optional)
%
% glob -
%
%   (limited) glob search string used to restrict which 
%   testscripts are run ("limited" in the sense that we only
%   recognize a string literal and a tilde-negated string literal. 
%   See example usage below). If not specified, all testscripts 
%   are run.
%
% haltonerror,useqsubdebug -
%
%   == true, testing will halt (crash) if any script errors
%   == false (default), testing will ignore error and simply run 
%      next script UNLESS useqsubdebug == true, then we instead 
%      drop into aaq_qsub_debug (default: false)
%
% Example Usage
%
%   aa_test;           - run all scripts, use default settings
%   aa_test([],true)   - run all scripts, halt on error
%
% example glob usage
%
%   aa_test('aatest_ds000114_fmri') - only run aatest_ds000114_fmri.m
%   aa_test('~ds002737') - run all scripts except those with 
%                          "ds002737" in the name (note leading tilde)
%
% test results (pass/fail) are logged to aa_test.log in working dir.
% If this file exists, you'll be asked before it is overwritten
%
% Expected usage is: 1) do 'aa_test' to run all scripts using
% defaults, 2) check the log, then 3) run aa_test('foo', true)
% or aa_test('foo',true,true) to re-run a script "foo" that failed
% in order to drop into the debugger.
%
% Notes
%
% test scripts are assumed to live in $AAHOME/developer/testscripts
% (old version assumed scripts lived in $AAHOME/examples). These
% scripts are designed for efficient PR testing and not necessarily
% helpful for learning aa. See $AAHOME/examples for the latter.
%
% Note testscripts load parameter file ~/.aa/aap_parameters_user.xml
% You should customize this file for how you want to do PR testing. 
% See comments in $AAHOME/developer/testscript_TEMPLATE.m for info.
%
% Revision History
%
% summer/2021 [MSJ] - newish (derived from aatest) 
%

if ~exist('haltonerror','var') || isempty(haltonerror)
    haltonerror = false;
end

if ~exist('useqsubdebug','var') || isempty(useqsubdebug)
    useqsubdebug = false;
end

if ~exist('glob','var') || isempty(glob)
    glob = '';
end

% parse glob negation tilde

globflag = 0;

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

aadir = fileparts(which('aa_ver5'));
testdir = fullfile(aadir,'developer/testscripts');

tests = dir(fullfile(testdir,'*.m'));
testnames = {tests.name};

savedir = pwd;
cd(testdir);

for findex = 1:numel(testnames)
    fname = testnames{findex};
    if globflag<0 && contains(fname,glob); continue; end
    if globflag>0 && ~contains(fname,glob); continue; end
    fprintf('\n\nRunning %s...\n', fname);
    runit(fullfile(testdir,fname), haltonerror, useqsubdebug, fid);
end

fprintf('------------------\ntests finished.\n')
cd(savedir);
fclose(fid);

end

% -------------------------------------------------------------------------
function runit(scriptpath, haltonerror, useqsubdebug, fid)
% -------------------------------------------------------------------------

if  haltonerror
    
    % run script normally;
    % function halts on aa error
    
    encapsulated(scriptpath);
    
    % if encapsulated returns, this script passed
    
    fprintf('\npass - %s\n',scriptpath);
    fprintf(fid,'\npass - %s\n',scriptpath);

else
    
    % run script inside a try-catch;
    % function flags any error and returns
    % optionally dropping into qsub_debug
    
    try
        
        encapsulated(scriptpath);
        fprintf('\npass - %s\n',scriptpath);
        fprintf(fid,'\npass - %s\n',scriptpath);
        
    catch err
        
        fprintf('FAIL - %s\n',scriptpath);
        fprintf(fid,'FAIL - %s\n',scriptpath);
        
        if useqsubdebug
            
            try
                aaq_qsub_debug;
            catch debugerr
                fprintf('aaq_qsub_debug. See err variable.\n');
                fprintf(fid,'aaq_qsub_debug. See err variable.\n');
            end
            
            % catch cases when aaq_qsub_debug quietly 
            % returns (e.g., failed to submitjob)
            
            keyboard;
            
        end
        
    end

end

end


% -------------------------------------------------------------------------
function encapsulated(scriptpath)
% -------------------------------------------------------------------------

% run script inside a function to prevent clear 
% etc from breaking outer variable scope

run(scriptpath)

end




