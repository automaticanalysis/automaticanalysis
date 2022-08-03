function [aap, resp] = aamod_meeg_sourcereconstruction(aap,task,subj)

resp='';
EST = aas_getstreams(aap,'output'); EST = setdiff(EST,{'peak'},'stable');

switch task
    case 'report'
        modulename = aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name;

        fnimg = fullfile(aas_getsubjpath(aap,subj),['diagnostic_' modulename '_models.jpg']);
        if exist(fnimg,'file')
            aap = aas_report_add(aap,subj,'<h4>Models</h4>');
            aap=aas_report_addimage(aap,subj,fnimg); 
        end
        
        fnimg = fullfile(aas_getsubjpath(aap,subj),['diagnostic_' modulename '_sensors.jpg']);
        if exist(fnimg,'file')
            aap = aas_report_add(aap,subj,'<h4>Sensors</h4>');
            aap=aas_report_addimage(aap,subj,fnimg); 
        end
        
        fnimg = fullfile(aas_getsubjpath(aap,subj),['diagnostic_' modulename '_leadfieldrank.jpg']);
        if exist(fnimg,'file')
            aap = aas_report_add(aap,subj,'<h4>Leadfield rank deficiency</h4>');
            aap=aas_report_addimage(aap,subj,fnimg); 
        end                
        
    case 'doit'
        [~, FT] = aas_cache_get(aap,'fieldtrip');
        FT.load;
        FT.addExternal('spm12');
        
        prepareOnly = strcmp(aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name,'aamod_meeg_preparesourcereconstruction');
        
        %% obtain inputs
        dat = load(aas_getfiles_bystream(aap,'subject',subj,'segmentation'));
        seg = dat.seg;
        dat = load(aas_getfiles_bystream(aap,'subject',subj,'headmodel'));
        headmodel = dat.headmodel;
        cfg = [];
        cfg.headshape = headmodel;
        cfg.numvertices = 1000;
        plotHeadmodel = ft_prepare_mesh(cfg);
        dat = load(aas_getfiles_bystream(aap,'subject',subj,'sourcemodel'));
        sourcemodel = dat.sourcemodel;
        hasSensors = aas_stream_has_contents(aap,'sensors');
        
        % band-average
        bndcfg = [];
        bndcfg.parameter = {'pow'};
        bndcfg.bandspec = aas_getsetting(aap,'bandspecification');

        inpstream = aas_getstreams(aap,'input'); inpstream = inpstream{end};
        if aas_stream_has_contents(aap,'subject',subj,inpstream)
            infnames = cellstr(aas_getfiles_bystream(aap,'subject',subj,inpstream));
        elseif aas_stream_has_contents(aap,'meeg_session',[subj 1],inpstream)
            infnames = cellstr(aas_getfiles_bystream(aap,'meeg_session',[subj 1],inpstream));
        else
            aas_log(aap,true,['Stream ' inpstream ' not found']);
        end
        switch spm_file(infnames{1},'ext')
            case 'mat'
                dat = load(infnames{1});
                f = fieldnames(dat);
                data = ft_struct2single(dat.(f{1}));
                outfnames = spm_file(infnames,'prefix','source_');
            case 'set'
                infnames = infnames(strcmp(spm_file(infnames,'ext'),'set'));
                FT.unload;
                [~, EL] = aas_cache_get(aap,'eeglab');
                EL.load;
                EEG = pop_loadset('filepath',spm_file(infnames{1},'path'),'filename',spm_file(infnames{1},'filename'));
                EL.unload;
                FT.reload;
                data = ft_struct2single(eeglab2fieldtripER(EEG,'reorient',1));
                outfnames = spm_file(infnames,'prefix','source_','ext','mat');
            otherwise
                aas_log(aap,true,'Unsupported file format')
        end
        
        if ~hasSensors
            elec = data.elec;
            
            % plot headmodel and sourcemodel
            f = figure; hold on;
            ft_plot_mesh(plotHeadmodel, 'maskstyle', 'opacity', 'facecolor', 'black', 'facealpha', 0.25, 'edgecolor', [0.5 0.5 0.5],   'edgeopacity', 0.5);
            ft_plot_mesh(sourcemodel, 'maskstyle', 'opacity', 'facecolor', 'black', 'facealpha', 0.25, 'edgecolor', 'red',   'edgeopacity', 0.5);
            view(135,45);
            set(f,'position',[0,0,720 720]);
            set(f,'PaperPositionMode','auto');
            print(f,'-noui',fullfile(aas_getsubjpath(aap,subj),['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_models']),'-djpeg','-r300');
            close(f);
            
            %% realign electrodes
            settings = aas_getsetting(aap,'realignelectrodes');
            
            % create target mesh
            % - create hull
            cfg = [];
            cfg.tissue = {settings.target};
            cfg.method = 'hexahedral';
            mesh = ft_prepare_mesh(cfg,seg);
            
            % - create low-res mesh
            cfg = [];
            cfg.method = 'headshape';
            cfg.headshape = mesh;
            cfg.numvertices = 1000;
            mesh = ft_prepare_mesh(cfg);
            
            mesh = ft_convert_units(mesh, elec.unit);
            
            % realign
            switch settings.method
                case 'fieldtrip'
                    cfg = [];
                    cfg.method = 'headshape';
                    cfg.warp = 'traditional';
                    cfg.headshape = mesh;
                    elec_realigned = ft_electroderealign(cfg,elec);
                case 'spherefit'
                    sphElec = pcfitsphere(pointCloud(elec.elecpos),0.1);
                    sphTarget = pcfitsphere(pointCloud(mesh.pos),2);
                    transP = [sphTarget.Center-sphElec.Center 0 0 0 sphTarget.Radius/sphElec.Radius*ones(1,3)];
                    elec_realigned = elec;
                    elecpos = ft_warp_apply(spm_matrix(transP), elec.elecpos, 'homogenous');
                    elec_realigned.elecpos = elecpos;
                    elec_realigned.chanpos = elecpos;
            end
            
            
            if settings.projecttotarget
                cfg = [];
                cfg.method = 'project';
                cfg.headshape = mesh;
                elec_final = ft_electroderealign(cfg,elec_realigned);
            else
                elec_final = elec_realigned;
            end
            
            % report
            aas_log(aap,false,sprintf('INFO: median distance between the electrodes and the mesh "%s" (original, realigned, final): %1.3f%s -> %1.3f%s -> %1.3f%s',aas_getsetting(aap,'realignelectrodes.target'),...
                ft_warp_error([],elec.elecpos, mesh),elec.unit,...
                ft_warp_error([],elec_realigned.elecpos, mesh),elec.unit,...
                ft_warp_error([],elec_final.elecpos, mesh),elec.unit...
                ));
            
            % plot headmodel and electrodes
            f = figure; hold on;
            ft_plot_mesh(plotHeadmodel, 'maskstyle', 'opacity', 'facecolor', 'black', 'facealpha', 0.25, 'edgecolor', [0.5 0.5 0.5],   'edgeopacity', 0.5);
            ft_plot_mesh(mesh.pos,'vertexcolor','skin')
            ft_plot_sens(elec_final,'style','red');
            view(135,45);
            set(f,'position',[0,0,720 720]);
            set(f,'PaperPositionMode','auto');
            print(f,'-noui',fullfile(aas_getsubjpath(aap,subj),['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_sensors']),'-djpeg','-r300');
            close(f);
            
            clear seg mesh elec elec_realigned sphElec sphTarget plotHeadmodel
            
            if prepareOnly
                sens = elec_final;
                outfname = fullfile(aas_getsubjpath(aap,subj),'sensors.mat');
                save(outfname,'sens');
                aap = aas_desc_outputs(aap,'subject',subj,'sensors',outfname);
            end
        else
            dat = load(aas_getfiles_bystream(aap,'subject',subj,'sensors'));
            elec_final = dat.sens;
        end
        
        %% source reconstruction
        settings = aas_getsetting(aap,['options.' aas_getsetting(aap,'method')]);
        hasLeadField = aas_stream_has_contents(aap,'leadfield');
        hasFilter = aas_stream_has_contents(aap,'sourcefilter');
        
        % compute the leadfield
        if ~hasLeadField && ~hasFilter
            if isfield(sourcemodel,'cfg') && isfield(sourcemodel.cfg,'urinside') % save template mask
                urinside = sourcemodel.cfg.urinside;
            else
                urinside = [];
            end
            cfg = [];
            cfg.sourcemodel = sourcemodel;
            cfg.headmodel   = headmodel;
            cfg.normalize = settings.normalize;
            data.elec = elec_final;
            sourcemodel = ft_prepare_leadfield(cfg, data);
            
            % QA leadfield
            rankLF = cellfun(@rank, sourcemodel.leadfield);
            ind = rankLF ~= mode(rankLF);
            if any(ind)
                aas_log(aap,false,sprintf('WARNING: low leadfield rank has been found in %d of %d location(s) -> they will be excluded',sum(ind),numel(ind)))
                
                % - create QA plot
                sourceQA = keepfields(sourcemodel,{'unit','cfg'});
                sourceQA.pos = sourcemodel.pos(ind,:);
                sourceQA.inside = sourcemodel.inside(ind);
                f = figure; hold on;
                ft_plot_mesh(sourcemodel, 'maskstyle', 'opacity', 'facecolor', [0.5 0.5 0.5], 'facealpha', 0.25, 'edgecolor', 'none',   'edgeopacity', 0.5);
                ft_plot_mesh(sourceQA, 'vertexcolor', 'red');
                view(135,45);
                set(f,'position',[0,0,720 720]);
                set(f,'PaperPositionMode','auto');
                print(f,'-noui',fullfile(aas_getsubjpath(aap,subj),['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_leadfieldrank']),'-djpeg','-r300');
                close(f);
                
                % - remove affected positions and tessels
                sourcemodel.pos(ind,:) = [];
                sourcemodel.inside(ind,:) = [];
                sourcemodel.leadfield(ind) = [];
                if isfield(sourcemodel,'tri')
                    [trind, trelem] = ndgrid(1:size(sourcemodel.tri,1),1:size(sourcemodel.tri,2)); trind = transpose(trind); trelem = transpose(trelem);
                    indtri = logical(sum(arrayfun(@(x,y) any(sourcemodel.tri(x,y)==find(ind)), trind, trelem)));
                    sourcemodel.tri = arrayfun(@(x,y) sourcemodel.tri(x,y)-sum(ind(1:sourcemodel.tri(x,y))), trind, trelem)';
                    sourcemodel.tri(indtri,:) = [];
                end
            end
            sourcemodel.included = reshape(~ind,[],1);
            if ~isempty(urinside)
                urind = find(urinside);
                urind(~sourcemodel.included) = [];
                sourcemodel.included = false(numel(urinside),1);
                sourcemodel.included(urind) = true;
            end
            
            clear rankLF sourceQA trind trelem indtri
            if prepareOnly
                outfname = fullfile(aas_getsubjpath(aap,subj),'leadfield.mat');
                save(outfname,'sourcemodel');
                aap = aas_desc_outputs(aap,'subject',subj,'leadfield',outfname);
            end
        end
        
        % run
        if hasLeadField
            dat = load(aas_getfiles_bystream(aap,'subject',subj,'leadfield'));
            sourcemodel = dat.sourcemodel;
        end
        if hasFilter
            dat = load(aas_getfiles_bystream(aap,'subject',subj,'sourcefilter'));
            sourcemodel.filter = dat.filter;
        end
        switch aas_getsetting(aap,'method')
            case 'eloreta'
                cfg             = [];
                % cfg.frequency   = ;
                cfg.method      = 'eloreta';
                cfg.sourcemodel = sourcemodel;
                cfg.headmodel   = headmodel;
                cfg.normalize = settings.normalize;
                cfg.eloreta.keepfilter = 'yes';
        end
        
        % Detect contrasts ansd recalculate
        sourceaap = aas_setcurrenttask(aap,aas_getsourcestage(aap,'aamod_meeg_timefrequencyanalysis'));
        models = aas_getsetting(sourceaap,'trialmodel');
        models = models(strcmp({models.subject},aas_getsubjname(aap,subj))).model;
        for i = 1:numel(infnames)
            m = models(arrayfun(@(m) endsWith(spm_file(infnames{i},'basename'),['_' m.name]), models));
            if prepareOnly
                data.elec = elec_final;
                timefreq     = ft_struct2single(ft_sourceanalysis(cfg, data));  % compute the source model
                filter = timefreq.avg.filter;
                outfname = fullfile(aas_getsubjpath(aap,subj),'filter.mat');
                save(outfname,'filter');
                aap = aas_desc_outputs(aap,'subject',subj,'sourcefilter',outfname);
                break;
            else
                dat = load(infnames{i}); f = fieldnames(dat); data = dat.(f{1});
                data.elec = elec_final;
                timefreq     = ft_struct2single(ft_sourceanalysis(cfg, data));  % compute the source model
            end
            timefreq.avg.dimord = strrep(data.dimord,'chan','pos');
            if numel(m.event.weights) == 1 % single
                if ~isempty(bndcfg.bandspec.band) % bands
                    ipf = [];
                    if aas_stream_has_contents(aap,'subject',subj,'ipf')
                        fnIPF = cellstr(aas_getfiles_bystream(aap,'subject',subj,'ipf'));
                        load(fnIPF{contains(spm_file(fnIPF,'basename'),m.event.names{1})},'ipf');
                    end
                    timeband = ft_average_bands(bndcfg,ipf,timefreq);
                end
            else  % contrast (assumed to be later than its primaries
                testmodel = repmat(rmfield(m,'name'),1,numel(m.event.names));
                for tm = 1:numel(testmodel)
                    testmodel(tm).event.names = m.event.names(tm);
                    testmodel(tm).event.weights = 1;
                end
                prmodels = arrayfun(@(tm) find(arrayfun(@(om) isequal(rmfield(om,'name'),tm),models)), testmodel, 'UniformOutput', false);
                if any(cellfun(@isempty, prmodels)), aas_log(aap,true,sprintf('One or more of the primary trialmodel(s) of %s not found', m.name)); end
                prmodels = cell2mat(prmodels);
                for estim = EST
                    for pri = 1:numel(prmodels)
                        dat = load(strrep(outfnames{pri},'timefreq',estim{1})); f = fieldnames(dat); dat = dat.(f{1});
                        if isfield(dat.avg,'dimord'), dat.avg = rmfield(dat.avg,'dimord'); end
                        prdata{pri} = ft_selectdata([],dat);
                        prdata{pri}.dimord = strrep(data.dimord,'chan','pos');
                    end
                    combinecfg.parameter = 'pow';
                    combinecfg.weights = m.event.weights;
                    combinecfg.contrast = aas_getsetting(sourceaap,'contrastoperation');
                    combinecfg.normalise = 'no';
                    combinedata = ft_combine(combinecfg,prdata{:});
                    switch estim{1}
                        case 'timefreq'
                            timefreq.avg.pow = combinedata.pow;
                        case 'timeband'
                            timeband.avg.pow = combinedata.pow;
                    end                    
                end
            end
            timefreq.avg = rmfield(timefreq.avg,intersect(fieldnames(timefreq.avg),{'filter','filterdimord'}));
            timefreq.cfg = []; % remove provenance to save space
            timefreq.cfg.included = sourcemodel.included;
            save(outfnames{i},'timefreq')
            if ~isempty(bndcfg.bandspec.band)
                timefreq = timeband;
                timefreq.avg.dimord = strrep_multi(data.dimord,{'chan' 'freq'},{'pos' 'band'});
                timefreq.avg = rmfield(timefreq.avg,intersect(fieldnames(timefreq.avg),{'filter','filterdimord'}));
                timefreq.cfg = []; % remove provenance to save space
                timefreq.cfg.included = sourcemodel.included;
                save(strrep(outfnames{i},'freq','band'),'timefreq')
            end
        end
        
        %% save outputs
        if ~prepareOnly
            aap = aas_desc_outputs(aap,'subject',subj,'timefreq',outfnames);
            if ~isempty(bndcfg.bandspec.band)
                aap = aas_desc_outputs(aap,'subject',subj,'timeband',strrep(outfnames,'freq','band'));
            end
        end
        
        FT.rmExternal('spm12');
        FT.unload;
    case 'checkrequirements'
        if ~aas_cache_get(aap,'fieldtrip'), aas_log(aap,true,'FieldTrip is not found'); end
        
        if strcmp(aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name,'aamod_meeg_sourcereconstruction')
            instream = aas_getstreams(aap,'input'); instream = instream{end};
            outstream = aas_getstreams(aap,'output',1); % assume single output
            instream = textscan(instream,'%s','delimiter','.'); instream = instream{1}{end};
            if ~strcmp(outstream,instream)
                aap = aas_renamestream(aap,aap.tasklist.currenttask.name,outstream,instream,'output');
                aas_log(aap,false,['INFO: ' aap.tasklist.currenttask.name ' output stream: ''' instream '''']);
                if ~isempty(aas_getsetting(aap,'bandspecification.band'))
                    aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'append',strrep(instream,'freq','band'),'output');
                    aas_log(aap,false,['INFO: ' aap.tasklist.currenttask.name ' output stream: ''' strrep(instream,'freq','band') '''']);
                end
            end
        end
end
end