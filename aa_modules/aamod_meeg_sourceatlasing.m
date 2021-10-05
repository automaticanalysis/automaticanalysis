function [aap, resp] = aamod_meeg_sourceatlasing(aap,task,varargin)

resp='';

switch task
    case 'report'

    case 'doit'
        [junk, FT] = aas_cache_get(aap,'fieldtrip');
        FT.load;
        FT.addExternal('spm12');
        
        instream = aas_getstreams(aap,'input'); instream = instream{end};
        inputfnames = cellstr(aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,cell2mat(varargin),instream));
        dat = load(inputfnames{1}); f = fieldnames(dat); source = dat.(f{1});
        
        %% Create atlas
        
        dat = load(aas_getfiles_bystream(aap,'subject',varargin{1},'sourceatlas')); sourceatlas = dat.sourceatlas;
        
        % resample
        cfg = [];
        cfg.interpmethod = 'nearest';
        cfg.parameter = 'aparc';
        source_atlas = ft_sourceinterpolate(cfg, sourceatlas, source);
        
        % remove unrepresented areas
        sel = unique(source_atlas.aparc); sel(sel==0) = [];
        source_atlas.aparclabel = source_atlas.aparclabel(sel);
        for s = 1:numel(sel), source_atlas.aparc(source_atlas.aparc==sel(s)) = s; end
        
        % create spatial information for 'sensors' N.B.: it is still suboptimal 
        for l = 1:numel(source_atlas.aparclabel)
            ind = source_atlas.aparc == l;
            labelpos(l,:) = mean(source_atlas.pos(ind,:),1);
        end
        elec.label = strrep(source_atlas.aparclabel,'_',' ');
        elec.chanpos = labelpos;
        elec.elecpos = labelpos;
        
        %% Run through inputs
        outputfnames = {};
        for ifn = 1:numel(inputfnames)
            dat = load(inputfnames{ifn}); f = fieldnames(dat); source = dat.(f{1});
            
            % Initialize label data
            labeldata = rmfield(source,intersect({'inside','pos','tri','method','avg'}, fieldnames(source)));
            labeldata.cfg = [];
            labeldata.elec = elec;
            labeldata.label = elec.label;
            
            if ~isfield(labeldata,'dimord')
                if isfield(source,'avg'), labeldata.dimord = source.avg.dimord; end
            end
            if isfield(labeldata,'dimord')
                labeldata.dimord = strrep(labeldata.dimord,'pos','label');
            else
                aas_log(aap,false,'WARNING: dimord not found');
            end
            
            % Atlasing
            for par = cellstr(aas_getsetting(aap,'parameter'))
                if ~isfield(source,par{1}) % look for parameter
                    if isfield(source,'avg'), source.(par{1}) = source.avg.(par{1}); end
                end
                if ~isfield(source,par{1})
                    aas_log(aap,false,sprintf('WARNING: parameters %s not found -> skipping',par{1}));
                    continue;
                end
                if iscell(source.(par{1})) %e.g. trials
                    parameter = {};
                    for t = 1:numel(source.(par{1}))
                        for l = 1:numel(source_atlas.aparclabel)
                            ind = source_atlas.aparc == l;
                            parameter{t}(l,:) = mean(source.(par{1}){t}(ind,:),1);
                        end
                    end
                else
                    parameter = [];
                    for l = 1:numel(source_atlas.aparclabel)
                        ind = source_atlas.aparc == l;
                        parameter(l,:) = mean(source.(par{1})(ind,:),1);
                    end
                end
                labeldata.(par{1}) = parameter;
            end
            outputfnames{ifn} = spm_file(inputfnames{ifn},'basename',strrep(spm_file(inputfnames{ifn},'basename'),'source','label'));
            save(outputfnames{ifn},'labeldata');
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