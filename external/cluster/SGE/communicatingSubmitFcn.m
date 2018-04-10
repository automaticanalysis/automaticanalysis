function communicatingSubmitFcn(cluster, job, props, varargin)
%COMMUNICATINGSUBMITFCN Submit a communicating MATLAB job to a SGE cluster
%
% Set your cluster's CommunicatingSubmitFcn to this function using the following
% command:
%     set(cluster, 'CommunicatingSubmitFcn', @communicatingSubmitFcn);
%
% See also parallel.cluster.generic.communicatingDecodeFcn.
%

% Copyright 2010-2015 The MathWorks, Inc.

% parse submission arguments
ind = find(cellfun(@(x) strcmp(x,'memory'),varargin));
if ind, req_memory = varargin{ind+1}; else req_memory = 2; end % default 2 GB
ind = find(cellfun(@(x) strcmp(x,'walltime'),varargin));
if ind, req_walltime = varargin{ind+1}; else req_walltime = 2; end % default 2 hours

% Store the current filename for the errors, warnings and dctSchedulerMessages
currFilename = mfilename;
if ~isa(cluster, 'parallel.Cluster')
    error('parallelexamples:GenericSGE:SubmitFcnError', ...
        'The function %s is for use with clusters created using the parcluster command.', currFilename)
end

decodeFunction = 'parallel.cluster.generic.communicatingDecodeFcn';

if ~cluster.HasSharedFilesystem
    error('parallelexamples:GenericSGE:SubmitFcnError', ...
        'The submit function %s is for use with shared filesystems.', currFilename)
end

if ~strcmpi(cluster.OperatingSystem, 'unix')
    error('parallelexamples:GenericSGE:SubmitFcnError', ...
        'The submit function %s only supports clusters with unix OS.', currFilename)
end

% The job specific environment variables
% Remove leading and trailing whitespace from the MATLAB arguments
matlabArguments = strtrim(props.MatlabArguments);
variables = {'MDCE_DECODE_FUNCTION', decodeFunction; ...
    'MDCE_STORAGE_CONSTRUCTOR', props.StorageConstructor; ...
    'MDCE_JOB_LOCATION', props.JobLocation; ...
    'MDCE_MATLAB_EXE', props.MatlabExecutable; ...
    'MDCE_MATLAB_ARGS', matlabArguments; ...
    'MDCE_DEBUG', 'true'; ...
    'MLM_WEB_LICENSE', props.UseMathworksHostedLicensing; ...
    'MLM_WEB_USER_CRED', props.UserToken; ...
    'MLM_WEB_ID', props.LicenseWebID; ...
    'MDCE_LICENSE_NUMBER', props.LicenseNumber; ...
    'MDCE_STORAGE_LOCATION', props.StorageLocation; ...
    'MDCE_CMR', cluster.ClusterMatlabRoot; ...
    'MDCE_TOTAL_TASKS', num2str(props.NumberOfTasks)};
% Set each environment variable to newValue if currentValue differs.
% We must do this particularly when newValue is an empty value,
% to be sure that we clear out old values from the environment.
for ii = 1:size(variables, 1)
    variableName = variables{ii,1};
    currentValue = getenv(variableName);
    newValue = variables{ii,2};
    if ~strcmp(currentValue, newValue)
        setenv(variableName, newValue);
    end
end
% Trim the environment variables of empty values.
nonEmptyValues = cellfun(@(x) ~isempty(strtrim(x)), variables(:,2));
variables = variables(nonEmptyValues, :);
% Which variables do we need to forward for the job?  Take all
% those that we have set.
variablesToForward = variables(:,1);

% Deduce the correct quote to use based on the OS of the current machine
if ispc
    quote = '"';
else
    quote = '''';
end


% The script name is communicatingJobWrapper.sh
scriptName = 'communicatingJobWrapper.sh';
% The wrapper script is in the same directory as this file
dirpart = fileparts(mfilename('fullpath'));
quotedScriptName = sprintf('%s%s%s', quote, fullfile(dirpart, scriptName), quote);

% Choose a file for the output. Please note that currently, JobStorageLocation refers
% to a directory on disk, but this may change in the future.
logFile = cluster.getLogLocation(job);
quotedLogFile = sprintf('%s%s%s', quote, logFile, quote);

jobName = sprintf('Job%d', job.ID);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CUSTOMIZATION MAY BE REQUIRED %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose a number of slots per node to use 
% You may wish to customize this section to match your cluster, 
% for example if you wish to limit the number of slots that 
% can be used for a single job.
% "-pe <parellel environment name>" - specify parallel environment
% "-q <queue name>" - specify queue
% "-l <resource list>" - specify resources
%   node=4 - request 4 nodes
%   mem=4gb - request 1GB RAM
%   walltime=1:00:00 - request 24h
procsPerNode = 1;
numberOfNodes = ceil(props.NumberOfTasks/procsPerNode);
% You may also wish to supply additional submission arguments to
% the qsub command here.
additionalSubmitArgs  = sprintf('-pe matlab %d', numberOfNodes);

% memory and walltime
additionalSubmitArgs = sprintf('%s -l h_rss=%dG -l h_cpu=%d:00:00',additionalSubmitArgs, req_memory, req_walltime);

dctSchedulerMessage(4, '%s: Requesting %d nodes with %d processors per node', currFilename, ...
    numberOfNodes, procsPerNode);
dctSchedulerMessage(5, '%s: Generating command for task %i', currFilename, ii);
commandToRun = getSubmitString(jobName, quotedLogFile, quotedScriptName, ...
    variablesToForward, additionalSubmitArgs);   

% Now ask the cluster to run the submission command
dctSchedulerMessage(4, '%s: Submitting job using command:\n\t%s', currFilename, commandToRun);
try
    % Make the shelled out call to run the command.
    [cmdFailed, cmdOut] = system(commandToRun);
catch err
    cmdFailed = true;
    cmdOut = err.message;
end
if cmdFailed
    error('parallelexamples:GenericSGE:SubmissionFailed', ...
        'Submit failed with the following message:\n%s', cmdOut);
end

dctSchedulerMessage(1, '%s: Job output will be written to: %s\nSubmission output: %s\n', currFilename, logFile, cmdOut);

jobIDs = extractJobId(cmdOut);
% jobIDs must be a cell array
if isempty(jobIDs)
    warning('parallelexamples:GenericSGE:FailedToParseSubmissionOutput', ...
        'Failed to parse the job identifier from the submission output: "%s"', ...
        cmdOut);
end
if ~iscell(jobIDs)
    jobIDs = {jobIDs};
end

% set the job ID on the job cluster data
cluster.setJobClusterData(job, struct('ClusterJobIDs', {jobIDs}));
