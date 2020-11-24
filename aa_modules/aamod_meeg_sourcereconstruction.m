function [aap, resp] = aamod_meeg_sourcereconstruction(aap,task,subj)

resp='';

switch task
    case 'report'

    case 'doit'
        [junk, FT] = aas_cache_get(aap,'fieldtrip');
        FT.load;
        FT.addExternal('spm12');
        
        %% obtain inputs
        mrifn = aas_getfiles_bystream(aap,'subject',subj,'structural');
        if strcmp(spm_file(mrifn,'ext'),'mat')
            dat = load(mrifn);
            mri = dat.mri;
        else
            aas_log(aap,true,'structural MUST be a Fieldtrip "mri"')
        end
        
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

        infnames = cellstr(aas_getfiles_bystream(aap,'subject',subj,'timefreq'));
        outfnames = spm_file(infnames,'prefix','source_');
        dat = load(infnames{1});
        elec = dat.timefreq.elec;
        
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
        
        %% source reconstruction
        settings = aas_getsetting(aap,['options.' aas_getsetting(aap,'method')]);
        
        % compute the leadfield
        cfg = [];
        cfg.sourcemodel = sourcemodel;
        cfg.headmodel   = headmodel;
        cfg.normalize = settings.normalize;
        dat.timefreq.elec = elec_final;
        sourcemodel = ft_prepare_leadfield(cfg, dat.timefreq);
                
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
        
        clear rankLF sourceQA trind trelem indtri
        
        % run
        switch aas_getsetting(aap,'method')
            case 'eloreta'
                cfg             = [];
                % cfg.frequency   = ;
                cfg.method      = 'eloreta';
                cfg.sourcemodel = sourcemodel;
                cfg.headmodel   = headmodel;
                cfg.normalize = settings.normalize;
        end
        
        for f = 1:numel(infnames)
            dat = load(infnames{f});
            dat.timefreq.elec = elec_final;
            timefreq     = ft_struct2single(ft_sourceanalysis(cfg, dat.timefreq));  % compute the source model
            timefreq.avg.dimord = strrep(dat.timefreq.dimord,'chan','pos');
            timefreq.cfg = []; % remove provenance to save space
            timefreq.cfg.included = reshape(~ind,[],1);
            save(outfnames{f},'timefreq')
        end

        %% save outputs        
        aap = aas_desc_outputs(aap,'subject',subj,'timefreq',outfnames);

        FT.rmExternal('spm12');
        FT.unload;
    case 'checkrequirements'
        if ~aas_cache_get(aap,'fieldtrip'), aas_log(aap,true,'FieldTrip is not found'); end
end
end