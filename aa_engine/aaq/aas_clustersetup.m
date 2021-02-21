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
pvpmod(varargin, {'mem', 'walltime'})

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
aas_log(aap,false,sprintf('INFO: pool %s is detected', clusterTypeString));
% - maxfilt module in tasklist?
isMaxFiltInTasklist = any(strcmp({aap.tasklist.main.module.name},'aamod_meg_maxfilt'));
% - neuromag specified?
isNeuromagSpec = ~isempty(aap.directory_conventions.neuromagdir);

switch class(obj.pool)
    case 'parallel.cluster.Slurm'
        obj.pool.SubmitArguments = sprintf('--mem=%dG --time=%d', mem, walltime*60);
        if isMaxFiltInTasklist && isNeuromagSpec
            obj.initialSubmitArguments = ' --constraint=maxfilter';
        end
        obj.pool.SubmitArguments = strcat(obj.pool.SubmitArguments, obj.initialSubmitArguments);
        if doSetNumWorkers
            aaparallel.numberofworkers = 1;
        end
        
    case 'parallel.cluster.Torque'
        obj.pool.SubmitArguments = sprintf('-l mem=%dGB, walltime=%d:00:00', mem, walltime);
        if isMaxFiltInTasklist && isNeuromagSpec
            % TODO: clarify whether/how NODESET should be eliminated
            obj.initialSubmitArguments = ' -W x=\"NODESET:ONEOF:FEATURES:MAXFILTER\"';
        end
        obj.pool.SubmitArguments = strcat(obj.pool.SubmitArguments, obj.initialSubmitArguments);
        if doSetNumWorkers
            aaparallel.numberofworkers = 1;
        end
        
    case 'parallel.cluster.LSF'
        obj.pool.SubmitArguments = sprintf(' -c %d -M %d -R "rusage[mem=%d:duration=%dh]"',...
            walltime*60, mem*1000, mem*1000, walltime);
        obj.pool.SubmitArguments = strcat(obj.initialSubmitArguments, obj.pool.SubmitArguments);
        if doSetNumWorkers
            aaparallel.numberofworkers = aap.options.aaparallel.numberofworkers;
        end
        
    case 'parallel.cluster.Generic'
        % TODO: starting with R2017b, AdditionalProperties is available, so
        % in reasonably recent Matlab versions there should be no need to
        % check for newGenericVersion
        if ~isprop(obj.pool.AdditionalProperties,'AdditionalSubmitArgs')
            aas_log(aap,false,'WARNING: Property "AdditionalSubmitArgs" not found.');
            aas_log(aap,false,'    "AdditionalSubmitArgs" must be listed within AdditionalProperties in the cluster profile in order to customise resource requirement and consequential queue selection.');
            aas_log(aap,false,'    Your jobs will be submitted to th default queue.');
        else
            obj.pool.AdditionalProperties.AdditionalSubmitArgs = sprintf('%s -l s_cpu=%d:00:00 -l s_rss=%dG', ...
                obj.initialSubmitArguments, walltime, mem);
        end
        if doSetNumWorkers
            aaparallel.numberofworkers = 1;
        end
        
    case 'parallel.cluster.Local'
        if doSetNumWorkers
            aaparallel.numberofworkers = aap.options.aaparallel.numberofworkers;
        end
end
