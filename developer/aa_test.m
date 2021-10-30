function aa_test(varargin)
%
% run /developer testscripts (new version -- replaces aatest)
%
% inputs (optional)
%
% 'glob', glob -
%   (limited) glob search string used to restrict which 
%   testscripts are run ("limited" in the sense that we only
%   recognize a string literal and a tilde-negated string literal. 
%   See example usage below). If not specified, all testscripts 
%   are run.
%
% 'deleteprevious', deleteprevious -
%   == true, testing will delete any previous tests to trigger re-execution
%      (default: true)
%
% 'haltonerror', haltonerror -
%   == true, testing will halt (crash) if any errors (default: false)
%
% 'wheretoprocess', wheretoprocess- 
%   == 'localsingle' (default), jobs within testing will run locally in a
%      serial fashion UNLESS wheretoprocess specifies differently, then
%      jobs instead will be submitted to the cluster as specified
%
% Example Usage
%
%   aa_test;                      - run all scripts, use default settings
%   aa_test('haltonerror',true)   - run all scripts, halt on error
%
% example glob usage
%
%   aa_test('glob','aatest_ds000114_fmri') - only run aatest_ds000114_fmri.m
%   aa_test('glob','~ds002737') - run all scripts except those with 
%                          "ds002737" in the name (note leading tilde)
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
% Notes
%
% test scripts are assumed to live in $AAHOME/developer/testscripts
% (old version assumed scripts lived in $AAHOME/examples). These
% scripts are designed for efficient PR testing and not necessarily
% helpful for learning aa. See $AAHOME/examples for the latter.
%
% Note testscripts load parameter file ~/.aa/aap_parameters_user.xml
% You should customize this file for how you want to do PR testing. 
% See comments in $AAHOME/developer/aatest_ds000114_TEMPLATE.m for info.
%
% Revision History
%
% summer/2021 [MSJ] - newish (derived from aatest) 
%

argParse = inputParser;
argParse.addParameter('glob','',@ischar);
argParse.addParameter('deleteprevious',true,@(x) islogical(x) || isnumeric(x));
argParse.addParameter('haltonerror',false,@(x) islogical(x) || isnumeric(x));
argParse.addParameter('wheretoprocess','localsingle',@ischar);
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
    runit(fname{1}, argParse.Results.deleteprevious, argParse.Results.haltonerror, argParse.Results.wheretoprocess, fid);
end

fprintf('------------------\ntests finished.\n')
cd(savedir);
fclose(fid);

end

% -------------------------------------------------------------------------
function runit(scriptname, deleteprevious, haltonerror, wheretoprocess, fid)
% -------------------------------------------------------------------------

func = str2func(strrep(scriptname,'.m',''));

if  haltonerror
    
    % run script normally;
    % function halts on aa error
    
    func(deleteprevious, wheretoprocess);
    
    % if encapsulated returns, this script passed
    
    fprintf('\npass - %s\n',scriptname);
    fprintf(fid,'\npass - %s\n',scriptname);

else
    
    % run script inside a try-catch;
    % function flags any error and returns
    % optionally dropping into qsub_debug
    
    try
        
        func(deleteprevious, wheretoprocess);
        fprintf('\npass - %s\n',scriptname);
        fprintf(fid,'\npass - %s\n',scriptname);
        
    catch err
        
        fprintf('FAIL - %s\n',scriptname);
        fprintf(fid,'FAIL - %s\n',scriptname);
        
        if strcmp(wheretoprocess,'localsingle')
            reporterror(scriptname,err);
        else
            
            try
                aaq_qsub_debug;
            catch err
                reporterror(scriptname,err);
            end
            
            % catch cases when aaq_qsub_debug quietly 
            % returns (e.g., failed to submitjob)
            
            keyboard;
            
        end
        
    end

end

end

function reporterror(scriptname,err)
msg = sprintf('%s had an error: %s\n',scriptname,err.message);
for e = 1:numel(err.stack)
    % Stop tracking to internal
    if strfind(err.stack(e).file,'distcomp'), break, end
    msg = [msg sprintf('<a href="matlab: opentoline(''%s'',%d)">in %s (line %d)</a>\n', ...
        err.stack(e).file, err.stack(e).line,...
        err.stack(e).file, err.stack(e).line)];
end
fprintf(msg)
end