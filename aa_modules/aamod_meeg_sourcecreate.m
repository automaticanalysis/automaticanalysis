function [aap, resp] = aamod_meeg_sourcecreate(aap,task,varargin)

resp='';

switch task
    case 'report'

    case 'doit'
        [junk, FT] = aas_cache_get(aap,'fieldtrip');
        FT.load;
        FT.addExternal('spm12');
        
        %% Obtain data
        clear data
        instream = aas_getstreams(aap,'input'); instream = instream{end};
        inputfnames = cellstr(aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,cell2mat(varargin),instream));
        switch spm_file(inputfnames{1},'ext')
            case 'mat'
                filetype = 'fieldtrip';
            case 'set'
                filetype = 'eeglab';
                inputfnames = inputfnames(strcmp(spm_file(inputfnames,'ext'),'set'));
            otherwise
                aas_log(aap,true,'Unsupported file format')
        end
        for seg = 1:numel(inputfnames)
            switch filetype
                case 'fieldtrip'
                    dat = load(inputfnames{seg});
                    data(seg) = ft_struct2single(dat.data);
                case 'eeglab'
                    FT.unload;
                    if seg == 1
                        [junk, EL] = aas_cache_get(aap,'eeglab');
                        EL.load;
                    else
                        EL.reload;
                    end
                    EEG = pop_loadset(inputfnames{seg});
                    if isempty(EEG.epoch)
                        aas_log(aap,false,sprintf('WARNING: segment # %d has no trial --> skipped',seg));
                        continue;
                    end
                    EL.unload;
                    FT.reload;
                    data(seg) = ft_struct2single(eeglab2fieldtripER(EEG,'reorient',1));
            end
        end
        data(cellfun(@isempty, {data.trial})) = []; % remove skipped segments
        inputfnames(cellfun(@isempty, {data.trial})) = []; 
        
        dat = load(aas_getfiles_bystream(aap,'subject',varargin{1},'sourcefilter')); filter = dat.filter;
        dat = load(aas_getfiles_bystream(aap,'subject',varargin{1},'leadfield')); source = keepfields(dat.sourcemodel,{'pos','tri','inside','included'});
        
        %% Run through inputs
        doParallel = true;
        aapoolprofile = strsplit(aap.directory_conventions.poolprofile,':'); poolprofile = aapoolprofile{1};
        if ~strcmp(aap.options.wheretoprocess,'qsub'), aas_log(aap,false,sprintf('WARNING: pool profile %s is not used via DCS/MPaS; therefore it may not work for parfor',poolprofile)); end
        try
            cluster = parcluster(poolprofile);
            if numel(aapoolprofile) > 1, cluster.ResourceTemplate = strjoin({aapoolprofile{2} cluster.ResourceTemplate}, ' '); end
            global aaworker;
            wdir = spm_file(tempname,'basename'); wdir = fullfile(aaworker.parmpath,wdir(1:8));
            aas_makedir(aap,wdir);
            cluster.JobStorageLocation = wdir;
            nWorkers = min([cluster.NumWorkers numel(data)]);
            pool = gcp('nocreate');
            if isempty(pool) || pool.NumWorkers < nWorkers, delete(pool); pool = parpool(cluster,nWorkers); end
        catch E
            aas_log(aap,false,['WARNING: ' poolprofile ' could not been initialised - ' E.message ' --> parallelisation is disabled']);
            doParallel = false;
        end
        
        if doParallel
            parfor i = 1:numel(data)
                outputfnames{i} = spm_file(inputfnames{i},'prefix','source_','ext','mat');
                create_sourcedata(data(i),cellstr(aas_getsetting(aap,'parameter')),filter,source,outputfnames{i});
            end
            delete(pool);
        else
            for i = 1:numel(data)
                outputfnames{i} = spm_file(inputfnames{i},'prefix','source_','ext','mat');
                create_sourcedata(data(i),cellstr(aas_getsetting(aap,'parameter')),filter,source,outputfnames{i});
            end
        end
        
        aap = aas_desc_outputs(aap,aap.tasklist.currenttask.domain,cell2mat(varargin),instream,outputfnames);
        
        FT.rmExternal('spm12');
        FT.unload;
    case 'checkrequirements'
        if ~aas_cache_get(aap,'fieldtrip'), aas_log(aap,true,'FieldTrip is not found'); end
        
        instream = aas_getstreams(aap,'input'); instream = instream{end};
        [stagename, index] = strtok_ptrn(aap.tasklist.currenttask.name,'_0');
        stageindex = sscanf(index,'_%05d');
        outstream = aap.tasksettings.(stagename)(stageindex).outputstreams.stream; % assume single output -> char
        instream = textscan(instream,'%s','delimiter','.'); instream = instream{1}{end};
        if ~strcmp(outstream,instream)
            aap = aas_renamestream(aap,aap.tasklist.currenttask.name,outstream,instream,'output');
            aas_log(aap,false,['INFO: ' aap.tasklist.currenttask.name ' output stream: ''' instream '''']);
        end
end
end
function create_sourcedata(data,parameters,filter,source,savepath)
sourcedata = rmfield(data,intersect(fieldnames(data),{'elec','label','hdr'}));
sourcedata.inside = source.inside;
sourcedata.pos = source.pos;
sourcedata.tri = source.tri;
sourcedata.included = source.included;
sourcedata.cfg = [];
for par = parameters
    fprintf('INFO: computing %s\n', par{1});
    fprintf('Calculating data for position: % 6d/% 6d',0,0);
    for pos = 1:numel(filter)
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b% 6d/% 6d',pos,numel(filter));
        
        if iscell(data.(par{1}))
            timeseries = zeros(3,0);
            for t = 1:numel(data.(par{1}))
                timeseries = [timeseries filter{pos} * data.(par{1}){t}];
            end
            u = svd(timeseries, 'econ');
            for t = 1:numel(data.(par{1}))
                parameter{t}(pos,:) = u(:,1)' * filter{pos} * data.(par{1}){t};
            end
        else
            timeseries = filter{pos} * data.(par{1});
            u = svd(timeseries, 'econ');
            parameter(pos,:) = u(:,1)' * timeseries;
        end
    end
    fprintf(' - Done.\n');
    sourcedata.(par{1}) = parameter;
end
save(savepath,'sourcedata');
end