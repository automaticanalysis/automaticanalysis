function obj = aas_clustersetup(obj, aap, doSetNumWorkers, NameValueArgs)
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

arguments
    obj
    aap struct
    doSetNumWorkers logical
    NameValueArgs.mem
    NameValueArgs.walltime
end

global aaparallel

% assign value to variables mem and walltime: if corresponding name-value
% pair was specified, use the corresponding value, otherwise pick values of
% globals
if isfield(NameValueArgs, 'mem')
    mem = NameValueArgs.mem;
else
    mem = aaparallel.memory;
end
if isfield(NameValueArgs, 'walltime')
    walltime = NameValueArgs.walltime;
else
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

% NOTE: original code version temporarily left in place for comparison

switch class(obj.pool)
    case 'parallel.cluster.Slurm'
        %         aas_log(obj.aap,false,'INFO: pool Slurm is detected');
        %         if isprop(obj.pool,'ResourceTemplate')
        %             obj.pool.ResourceTemplate = sprintf('--ntasks=^N^ --cpus-per-task=^T^ --mem=%dG --time=%d', aaparallel.memory,aaparallel.walltime*60);
        %         else
        %             obj.pool.SubmitArguments = sprintf('--mem=%dG --time=%d', aaparallel.memory,aaparallel.walltime*60);
        %         end
        %         if any(strcmp({aap.tasklist.main.module.name},'aamod_meg_maxfilt')) && ... % maxfilt module detected
        %                 ~isempty(aap.directory_conventions.neuromagdir) % neuromag specified
        %             obj.initialSubmitArguments = ' --constraint=maxfilter';
        %         end
        %         obj.pool.SubmitArguments = strcat(obj.pool.SubmitArguments,obj.initialSubmitArguments);
        %         aaparallel.numberofworkers = 1;
        obj.pool.SubmitArguments = sprintf('--mem=%dG --time=%d', mem, walltime*60);
        if isMaxFiltInTasklist && isNeuromagSpec
            obj.initialSubmitArguments = ' --constraint=maxfilter';
        end
        obj.pool.SubmitArguments = strcat(obj.pool.SubmitArguments, obj.initialSubmitArguments);
        if doSetNumWorkers
            aaparallel.numberofworkers = 1;
        end
        
    case 'parallel.cluster.Torque'
        %         aas_log(obj.aap,false,'INFO: pool Torque is detected');
        %         obj.pool.ResourceTemplate = sprintf('-l nodes=^N^,mem=%dGB,walltime=%d:00:00', aaparallel.memory,aaparallel.walltime);
        %         if any(strcmp({aap.tasklist.main.module.name},'aamod_meg_maxfilt')) && ... % maxfilt module detected
        %                 ~isempty(aap.directory_conventions.neuromagdir) % neuromag specified
        %             obj.initialSubmitArguments = ' -W x=\"NODESET:ONEOF:FEATURES:MAXFILTER\"';
        %         end
        %         obj.pool.SubmitArguments = strcat(obj.pool.SubmitArguments,obj.initialSubmitArguments);
        %         aaparallel.numberofworkers = 1;
        obj.pool.SubmitArguments = sprintf('mem=%dGB, walltime=%d:00:00', mem, walltime);
        if isMaxFiltInTasklist && isNeuromagSpec
            % TODO: clarify whether/how NODESET should be eliminated
            obj.initialSubmitArguments = ' -W x=\"NODESET:ONEOF:FEATURES:MAXFILTER\"';
        end
        obj.pool.SubmitArguments = strcat(obj.pool.SubmitArguments, obj.initialSubmitArguments);
        if doSetNumWorkers
            aaparallel.numberofworkers = 1;
        end
        
    case 'parallel.cluster.LSF'
        %         aas_log(obj.aap,false,'INFO: pool LSF is detected');
        %         obj.pool.SubmitArguments = sprintf(' -c %d -M %d -R "rusage[mem=%d:duration=%dh]"',aaparallel.walltime*60, aaparallel.memory*1000,aaparallel.memory*1000,aaparallel.walltime);
        %         obj.pool.SubmitArguments = strcat(obj.initialSubmitArguments,obj.pool.SubmitArguments);
        %         aaparallel.numberofworkers = aap.options.aaparallel.numberofworkers;
        obj.pool.SubmitArguments = sprintf(' -c %d -M %d -R "rusage[mem=%d:duration=%dh]"',...
            walltime*60, mem*1000, mem*1000, walltime);
        obj.pool.SubmitArguments = strcat(obj.initialSubmitArguments, obj.pool.SubmitArguments);
        if doSetNumWorkers
            aaparallel.numberofworkers = aap.options.aaparallel.numberofworkers;
        end
        
    case 'parallel.cluster.Generic'
        %         aas_log(obj.aap,false,'INFO: Generic engine is detected');
        %         obj.newGenericVersion = isempty(obj.pool.IndependentSubmitFcn);
        %         if obj.newGenericVersion
        %             if ~isprop(obj.pool.AdditionalProperties,'AdditionalSubmitArgs')
        %                 aas_log(obj.aap,false,'WARNING: Propertiy "AdditionalSubmitArgs" not found.');
        %                 aas_log(obj.aap,false,'    "AdditionalSubmitArgs" must be listed within AdditionalProperties in the cluster profile in order to customise resource requirement and consequential queue selection.');
        %                 aas_log(obj.aap,false,'    Your jobs will be submitted to th default queue.');
        %             else
        %                 obj.pool.AdditionalProperties.AdditionalSubmitArgs = sprintf('%s -l s_cpu=%d:00:00 -l s_rss=%dG',obj.initialSubmitArguments,aaparallel.walltime,aaparallel.memory);
        %             end
        %         else
        %             obj.pool.IndependentSubmitFcn = obj.SetArg(obj.pool.IndependentSubmitFcn,'walltime',aaparallel.walltime);
        %             obj.pool.IndependentSubmitFcn = obj.SetArg(obj.pool.IndependentSubmitFcn,'memory',aaparallel.memory);
        %         end
        %         aaparallel.numberofworkers = 1;
        
        
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
        %         aas_log(obj.aap,false,'INFO: Local engine is detected');
        %         aaparallel.numberofworkers = aap.options.aaparallel.numberofworkers;
        if doSetNumWorkers
            aaparallel.numberofworkers = aap.options.aaparallel.numberofworkers;
        end
end
