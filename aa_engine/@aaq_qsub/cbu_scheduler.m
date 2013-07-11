function scheduler=cbu_scheduler(varargin)
% Configure a scheduler object to use with cbu queueing system
% User can choose from several pre-specified configurations, or specify
% their own, custom configuration.
%
% Input Arguments:
%
% scheduler name (optional; string)
%   Name of the scheduler configuration to use.
%   If no name is specified, the function returns a default configuration.
%
%   Available options:
%       - 'basic-compute', 'compute', 'basic'
%               compute queue, 12 workers, 4Gb RAM per worker [default config]
%
%       - 'large-compute'
%               compute queue, 24 workers, 12Gb RAM per worker
%
%       - 'basic-gpu', 'gpu'
%               gpu queue, 12 workers, 4Gb RAM
%           
%       - 'custom'
%               maunally specify scheduler parameters (see below)
%
% scheduler parameters (optional; 1xn cell array)
%   This argument is only used with 'custom' config.
% 	Values for scheduler parameters specified in the following order:
%
%       queue name (compute, gpu, super)    [default = 'compute']
%       n workers                           [default = 12]
%       memory (Gb)                         
%           [maximum physical memory required by ALL of the workers on each 
%            physical host. Default = 4Gb]
%       job data directory                  
%           [location matlab uses to store information for distributed 
%            compute jobs. Default = /imaging/<user name>/.cbu-cluster/matlab-jobs]
%       matlab worker path
%           [which version of matlab should be used for workers. Default =
%           matlabroot - i.e. workers will use the same version of matlab
%           as is used to submit the job]
%
%   Any of these parameters can be left unspecified, in which case they will
%   be set to the default value. Parameters must be specified in the correct
%   order, however.
%
%
% queue submission string (optional; string)
%   This argument is only used with 'custom' config.
%   It provides a way to enter a completely custom queue submission string,
%   specifying submission arguments in the form used by the qsub command. 
%   There is no default value. 
%   If a queue sub string is specified, it will override any values entered
%   for queue name and memory in the scheduler parameters array (i.e. those
%   values will be ignored)

%
% Usage:
% e.g.  
%       scheduler = cbu_scheduler;
%       scheduler=cbu_scheduler('basic');
%       scheduler=cbu_scheduler('custom',{'gpu', 16, 64});
%       scheduler=cbu_scheduler('custom',{[],[], 64});
%       scheduler=cbu_scheduler('custom',{[],24,[],[],'/hpc-software/matlab/r2009a'},'-q gpu -l mem=96gb');
%
%       cbu_qsub(my_jobs,cbu_scheduler('gpu'));

if nargin <1
    varargin{1}='compute';
end


switch varargin{1}
    
    case {'basic-compute', 'compute', 'basic'}
         scheduler=format_sched({'compute',12,4});
         
    case {'large-compute'}
         scheduler=format_sched({'compute',24,12});
         
    case {'basic-gpu','gpu'}
         scheduler=format_sched({'gpu',12,4});
         
    case 'custom'
         if nargin<2
                error('cbu_scheduler:Missing input argument - you have specified a custom configuration but not provided the custom parameters');
         end
         scheduler=format_sched(varargin{2:end});
        
    otherwise
        error(['cbu_scheduler: Scheduler name ' varargin{1} ' not recognised']);
        
end


function sched=format_sched(varargin)

% get scheduler
if isempty(which('parcluster'))
    % pre R2012a, parcluster is not available, use the old findResource
    % method
    sched=findResource('scheduler','type','torque'); %#ok<REMFF1> % #ok<MSNU>  %#ok Mlint
    
    % the different types of scheduler objects have different fieldnames...
    fieldnames.has_shared='HasSharedFilesystem';
    fieldnames.sub_arg='SubmitArguments';
    fieldnames.workers='ClusterSize';
    fieldnames.job_loc='DataLocation';
    fieldnames.matlab='ClusterMatlabRoot';
    
    job_dir = 'pre_r2012a';
else
    % new method R2012a onwards
    
    profiles=parallel.clusterProfiles;
    profiles=strfind(profiles,'CBU_Cluster');
    if isempty(vertcat(profiles{:}))
        parallel.importProfile('/hpc-software/matlab/cbu/CBU_Cluster.settings');
    end
    
    sched=parcluster('CBU_Cluster');
    
    % fieldnames for the older scheduler object
    fieldnames.has_shared='HasSharedFilesystem';
    fieldnames.sub_arg='SubmitArguments';
    fieldnames.workers='NumWorkers';
    fieldnames.job_loc='JobStorageLocation';
    fieldnames.matlab='ClusterMatlabRoot';
    
    job_dir = 'r2012a_onwards';
end


% set up scheduler parameters

% default parameters
def.queue='compute';
def.nworkers=12;
def.memory=4;
def.job_dir=get_job_dir(job_dir); 	% location matlab uses to store information for
                                    % distributed compute jobs. Default = 
                                    % /imaging/<user name>/.cbu-cluster/matlab-jobs
                            
def.matlab_root=matlabroot; % location of matlab install on worker nodes. 
                            % Default = matlabroot - i.e. the root directory 
                            % of the matlab version used to run the
                            % submission script. Could allow people to
                            % specify a different matlab version be used on
                            % worker nodes, but how often would this be
                            % used? In reality, you'd just submit using a
                            % different version of matlab.




% configure parameters
sched.(fieldnames.has_shared)=logical(1); %#ok<LOGL>



% queue submission arguments:
if nargin>1 && ~isempty(varargin{2})
    % if q sub string is specified, use that:
    % ## TO DO - error checking?
    sched.(fieldnames.sub_arg)=varargin{2};
else
    % build qsub string from individual arguments

    % queue name:
    if size(varargin{1},2)<1 || isempty(varargin{1}{1})
        queue = def.queue;
    else
        if ~strcmpi(varargin{1}{1},'compute') && ~strcmpi(varargin{1}{1},'gpu') && ~strcmpi(varargin{1}{1},'super')
            error(['cbu_scheduler: Specified queue name ' varargin{1}{1} ' is not recognised'])
        else
            queue = varargin{1}{1};
        end
    end

    % max memory required for each worker that will be launched on
    % each physical host.
    if size(varargin{1},2)<3 || isempty(varargin{1}{3})
        maxmemory = def.memory;
    else
       if isnumeric(varargin{1}{3})
           if varargin{1}{3} > 144
               error(['cbu_scheduler: The maximum memory value specified (' varargin{1}{3} 'Gb) is too large. There are currently no machines available with this amount of memory']);
           elseif varargin{1}{3} > 96
               disp(['cbu_scheduler: The maximum memory value specified (' varargin{1}{3} 'Gb) is very large. There is currently only one machine available with this amount of memory']);
           elseif varargin{1}{3} > 16
               disp(['cbu_scheduler: You have reserved a large amount of memory (' varargin{1}{3} 'Gb)' ...
                   ' for each worker node. This will mean you can only run a limited number of workers ' ...
                   'on each physical host, and may slow down execution of your job ']);
           end
           maxmemory =varargin{1}{3};
        else
            error('cbu_scheduler: The maximum memory value is not specified correctly (value given must be numeric)')
        end
    end
    
    % ## TO DO - other qsub options?
    
    % put them together:
    sched.(fieldnames.sub_arg) = ['-q ' queue ' -l mem=' num2str(maxmemory) 'gb'];
end



% n workers
if size(varargin{1},2)<2 || isempty(varargin{1}{2})
    sched.(fieldnames.workers) = def.nworkers;
else
    % limit max workers to 32
    if isnumeric(varargin{1}{2})
        sched.(fieldnames.workers) = min(varargin{1}{2},32);
    else
        error('cbu_scheduler: The number of workers is not specified correctly (value given must be numeric)')
    end
end



% Scheduler job directory
if size(varargin{1},2)<4 || isempty(varargin{1}{4})
    sched.(fieldnames.job_loc)=def.job_dir;
else
    if exist(varargin{1}{4},'dir')~=7
        disp(['cbu_scheduler: Specified job output directory ' varargin{1}{4} ' does not exist. This directory will be created.']);
        mkdir(varargin{1}{4});
    end
    sched.(fieldnames.job_loc)=varargin{1}{4};
end

% Worker node matlab root
if size(varargin{1},2)<5 || isempty(varargin{1}{5})
    sched.(fieldnames.matlab)=def.matlab_root;
else
    if exist(varargin{1}{5},'dir')~=7
        error(['cbu_scheduler: Specified Matlab root directory ' varargin{1}{5} ' does not exist']);
    end
    sched.(fieldnames.matlab)=varargin{1}{5};
end

    

function job_dir=get_job_dir(job_dir)
% get default matlab DCT data location
%/imaging/<user name>/.cbu-cluster/matlab-jobs

[rtn, whoami]=system('whoami');  %#ok<ASGLU> % #ok<MSNU> %#ok Mlint
% v2009b and later: [~, whoami]=system('whoami');
job_dir=['/imaging/' strtrim(whoami) '/.cbu-cluster/matlab-jobs_' job_dir];
if exist(job_dir,'dir')~=7
    mkdir(job_dir);
end

