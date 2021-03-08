function obj = aas_clustersetup(obj, aap, doSetNumWorkers, varargin)
% AAS_CLUSTERSETUP set up cluster object.
% 
% aas_clustersetup(obj, aap, false) sets up cluster object obj.pool (any of
%   Slurm | Torque | LSF | Generic | Local), using parameter values in 
%   global struct aaparallel.
% aas_clustersetup(obj, aap, true) sets up cluster object obj.pool in the
%   same manner, and additionally alters aaparallel.numberofworkers 
%   depending on the cluster type.
% aas_clustersetup(..., 'mem', 4) overrides aaparallel.memory, using the
%   specified value in GB.
% aas_clustersetup(..., 'walltime', 36) overrides aaparallel.walltime,
%   using the specified value in hours.

global aaparallel

mem = [];
walltime = [];
% deal with name-value input pairs
pvpmod(varargin, ["mem", "walltime"])

% assign value to variables mem and walltime: if corresponding name-value
% pair was specified, use the corresponding value, otherwise pick values of
% globals
if isempty(mem)
    mem = aaparallel.memory;
end
if isempty(walltime)
    walltime = aaparallel.walltime;
end

% common preparatory work:
% - log message 
clusterTypeString = erase(string(class(obj.pool)), "parallel.cluster.");
aas_log(aap, false, compose("INFO: pool %s is detected", clusterTypeString));

% - maxfilt module in tasklist?
isMaxFiltInTasklist = any(strcmp({aap.tasklist.main.module.name},"aamod_meg_maxfilt"));
% - neuromag specified?
isNeuromagSpec = ~isempty(aap.directory_conventions.neuromagdir);

switch class(obj.pool)
    % new: Matlab Job Scheduler
    case 'parallel.cluster.MJS'
        if doSetNumWorkers
            aaparallel.numberofworkers = aap.options.aaparallel.numberofworkers;
        end
    
    case 'parallel.cluster.Slurm'
        obj.pool.SubmitArguments = compose("--mem=%dG --time=%d", mem, walltime*60);
        if isMaxFiltInTasklist && isNeuromagSpec
            obj.initialSubmitArguments = " --constraint=maxfilter";
        end
        obj.pool.SubmitArguments = convertStringsToChars(obj.pool.SubmitArguments +...
            obj.initialSubmitArguments);
        if doSetNumWorkers
            aaparallel.numberofworkers = 1;
        end
        
    case 'parallel.cluster.Torque'
        obj.pool.SubmitArguments = compose("-l mem=%dGB, walltime=%d:00:00", mem, walltime);
        if isMaxFiltInTasklist && isNeuromagSpec
            % TODO: clarify whether/how NODESET should be eliminated
            obj.initialSubmitArguments = " -W x=\""NODESET:ONEOF:FEATURES:MAXFILTER\""";
        end
        obj.pool.SubmitArguments = convertStringsToChars(obj.pool.SubmitArguments + ...
            obj.initialSubmitArguments);
        if doSetNumWorkers
            aaparallel.numberofworkers = 1;
        end
        
    case 'parallel.cluster.LSF'
        obj.pool.SubmitArguments = compose(" -c %d -M %d -R ""rusage[mem=%d:duration=%dh]""",...
            walltime*60, mem*1000, mem*1000, walltime);
        obj.pool.SubmitArguments = convertStringsToChars(obj.initialSubmitArguments + ...
            obj.pool.SubmitArguments);
        if doSetNumWorkers
            aaparallel.numberofworkers = aap.options.aaparallel.numberofworkers;
        end
        
    case 'parallel.cluster.Generic'
        obj.pool.AdditionalProperties.AdditionalSubmitArgs = convertStringsToChars(...
            string(obj.initialSubmitArguments) + ...
            compose(" -l s_cpu=%d:00:00 -l s_rss=%dG", walltime, mem));
        if doSetNumWorkers
            aaparallel.numberofworkers = 1;
        end
        
    case 'parallel.cluster.Local'
        if doSetNumWorkers
            aaparallel.numberofworkers = aap.options.aaparallel.numberofworkers;
        end
end
