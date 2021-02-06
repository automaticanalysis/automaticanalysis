function obj = aas_clustersetup(obj, aap)
% AAS_CLUSTERSETUP set up cluster object.
% 
% aas_clustersetup(obj, aap) sets up cluster object obj.pool (any of Slurm
%   | Torque | LSF | Generic | Local).


global aaparallel

poolTypeString = erase(string(class(obj.pool)), "parallel.cluster.");
aas_log(aap,false,sprintf('INFO: pool %s is detected', poolTypeString));

switch class(obj.pool)
    case 'parallel.cluster.Slurm'
        if isprop(obj.pool,'ResourceTemplate')
            obj.pool.ResourceTemplate = sprintf('--ntasks=^N^ --cpus-per-task=^T^ --mem=%dG --time=%d', aaparallel.memory,aaparallel.walltime*60);
        else
            obj.pool.SubmitArguments = sprintf('--mem=%dG --time=%d', aaparallel.memory,aaparallel.walltime*60);
        end
        if any(strcmp({aap.tasklist.main.module.name},'aamod_meg_maxfilt')) && ... % maxfilt module detected
                ~isempty(aap.directory_conventions.neuromagdir) % neuromag specified
            obj.initialSubmitArguments = ' --constraint=maxfilter';
        end
        obj.pool.SubmitArguments = strcat(obj.pool.SubmitArguments,obj.initialSubmitArguments);
        aaparallel.numberofworkers = 1;
    case 'parallel.cluster.Torque'
        obj.pool.ResourceTemplate = sprintf('-l nodes=^N^,mem=%dGB,walltime=%d:00:00', aaparallel.memory,aaparallel.walltime);
        if any(strcmp({aap.tasklist.main.module.name},'aamod_meg_maxfilt')) && ... % maxfilt module detected
                ~isempty(aap.directory_conventions.neuromagdir) % neuromag specified
            obj.initialSubmitArguments = ' -W x=\"NODESET:ONEOF:FEATURES:MAXFILTER\"';
        end
        obj.pool.SubmitArguments = strcat(obj.pool.SubmitArguments,obj.initialSubmitArguments);
        aaparallel.numberofworkers = 1;
    case 'parallel.cluster.LSF'
        obj.pool.SubmitArguments = sprintf(' -c %d -M %d -R "rusage[mem=%d:duration=%dh]"',aaparallel.walltime*60, aaparallel.memory*1000,aaparallel.memory*1000,aaparallel.walltime);
        obj.pool.SubmitArguments = strcat(obj.initialSubmitArguments,obj.pool.SubmitArguments);
        aaparallel.numberofworkers = aap.options.aaparallel.numberofworkers;
    case 'parallel.cluster.Generic'
        obj.newGenericVersion = isempty(obj.pool.IndependentSubmitFcn);
        if obj.newGenericVersion
            if ~isprop(obj.pool.AdditionalProperties,'AdditionalSubmitArgs')
                aas_log(aap,false,'WARNING: Propertiy "AdditionalSubmitArgs" not found.');
                aas_log(aap,false,'    "AdditionalSubmitArgs" must be listed within AdditionalProperties in the cluster profile in order to customise resource requirement and consequential queue selection.');
                aas_log(aap,false,'    Your jobs will be submitted to th default queue.');
            else
                obj.pool.AdditionalProperties.AdditionalSubmitArgs = sprintf('%s -l s_cpu=%d:00:00 -l s_rss=%dG',obj.initialSubmitArguments,aaparallel.walltime,aaparallel.memory);
            end
        else
            obj.pool.IndependentSubmitFcn = obj.SetArg(obj.pool.IndependentSubmitFcn,'walltime',aaparallel.walltime);
            obj.pool.IndependentSubmitFcn = obj.SetArg(obj.pool.IndependentSubmitFcn,'memory',aaparallel.memory);
        end
        aaparallel.numberofworkers = 1;
    case 'parallel.cluster.Local'
        aaparallel.numberofworkers = aap.options.aaparallel.numberofworkers;
end
