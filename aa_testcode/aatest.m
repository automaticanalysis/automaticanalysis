% test AA example user scripts. If testtype is undefined, we use the examples in the
% root of [aadir]/examples. If testtype is specified, it refers to a sub-directory, e.g.
% [aadir]/examples/cbu. By default we attempt to debug any crashes using aaq_qsub_debug.
% If you set the skipdebug flag to true we attempt to carry on instead.
%
% TODO: 
% * handle clean-up of old runs. At the moment it's up to the user to ensure all old
% results are removed before running.
% * support cluster submission of user scripts instead of running tests in serial
%
% aatest(testtype,[skipdebug=0])
function aatest(testtype,skipdebug)

if ~exist('testtype','var') || isempty(testtype)
    testtype = '';
end

if ~exist('skipdebug','var') || isempty(skipdebug)
    skipdebug = false;
end

aadir = fileparts(fileparts(mfilename('fullpath')));

testdir = fullfile(aadir,'examples',testtype);
tests = dir(fullfile(testdir,'aa_user_*.m'));
testnames = {tests.name};
% detect 'connect' user scripts - these need to run in a second pass since they depend
% on some previous results being generated
connectind = cell_index(testnames,'_connect');
% cell_index should probably return empty when nothing is found, but it doesn't
connectind(connectind==0) = [];
secondpass = testnames(connectind);
firstpass = setdiff(testnames,secondpass);

orgdir = pwd;
cd(testdir);
for thistest = firstpass(:)'
    tstr = thistest{1};
    fprintf('running test %s...\n',tstr);
    runit(fullfile(testdir,tstr),skipdebug);
end
for thistest = secondpass(:)'
    tstr = thistest{1};
    fprintf('running test %s...\n',tstr);
    runit(fullfile(testdir,tstr),skipdebug);
end
fprintf('------------------\ntests finished.\n')
cd(orgdir);

function runit(scriptpath,skipdebug)

try
    encapsulated(scriptpath);
    fprintf('pass - %s\n',scriptpath);
catch err
    fprintf('FAIL - %s\n',scriptpath);
    if ~skipdebug
        try
            aaq_qsub_debug;
        catch debugerr
            fprintf('DEBUG ERROR - error before aa_doprocessing? see err variable.\n');
        end
        % catch cases when aaq_qsub_debug quietly returns (e.g., failed to submit
        % job)
        keyboard;
    end
end

function encapsulated(scriptpath)
% run script in function to prevent clear etc from breaking outer variable scope.
run(scriptpath)
