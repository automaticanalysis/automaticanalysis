function [aap, resp] = aamod_meeg_prepareheadmodel(aap,task,subj)

resp='';

switch task
    case 'report'

    case 'doit'
        [junk, FT] = aas_cache_get(aap,'fieldtrip');
        FT.load;
        FT.addExternal('spm12');
        
        %% Obtain structural
        mrifn = aas_getfiles_bystream(aap,'subject',subj,'structural');
        if strcmp(spm_file(mrifn,'ext'),'mat')
            dat = load(mrifn);
            mri = dat.mri;
        else
            aas_log(aap,true,'structural MUST be a Fieldtrip "mri"')
        end
        
        %% Prepare headmodel
        settings = aas_getsetting(aap,['options.' aas_getsetting(aap,'method')]);
        
        inpstr = aas_getstreams(aap,'input');
        
        tpms = {'gray', 'white', 'csf', 'bone', 'softtissue', 'air'};
        if numel(inpstr) == 8 && all(cellfun(@(x) aas_stream_has_contents(aap,x),inpstr(2:8)))
            % import segmentation
            aas_log(aap,false,'INFO: suitable segmentation (6-class) found -> importing segmentation')
            
            seg = rmfield(mri,{'anatomy','inside'});
            seg.dim = double(seg.dim);
            targfname = fullfile(aas_getsubjpath(aap,subj),'target.nii');
            ft_write_mri(targfname,mri.anatomy,'transform',mri.transform,'dataformat','nifti')
            
            inpfname = cellfun(@(x) aas_getfiles_bystream(aap,'subject',subj,x),inpstr(2:8),'UniformOutput',false);
            structfname = inpfname{1};
            
            flags.estimate.cost_fun = 'nmi';
            flags.write.which = [1 0];
            flags.write.interp = 4;
            flags.write.mask = 1;
            flags.write.prefix = 'r';
            
            x = spm_coreg(spm_vol(targfname), spm_vol(structfname), flags.estimate);
                        
            for s = 2:7
                spm_get_space(inpfname{s}, spm_matrix(x)\spm_get_space(inpfname{s}));
                spm_reslice(char(targfname,inpfname{s}),flags.write);
                fname = spm_file(inpfname{s},'prefix',flags.write.prefix);
                seg.(tpms{s-1}) = spm_read_vols(spm_vol(fname));
                delete(fname);
            end
            
            delete(targfname);
        else
            aas_log(aap,false,'INFO: no suitable segmentation (6-class) found -> segmenting with FieldTrip')
            % create segmentation
            cfg          = [];
            cfg.spmmethod = 'new';
            cfg.opts.biasreg = 0.001;
            cfg.opts.samp = 1;
            cfg.output = {'tpm'};
            seg = ft_volumesegment(cfg, mri);
        end
        
        % - correct for padding        
        for f = tpms
            seg.(f{1})(isnan(mri.anatomy) | (mri.anatomy==0)) = 0;
        end
        
        % - smooth segmentation
        if aas_getsetting(aap,'segmentation.smoothing')
            for f = tpms
                Y = seg.(f{1}) + 0;
                spm_smooth(Y, Y, aas_getsetting(aap,'segmentation.smoothing'));
                seg.(f{1}) = Y;
            end
        end

        % - threshold segmentation
        segAll = zeros(seg.dim);
        scalp = seg.softtissue; % save for later
        switch aas_getsetting(aap,'segmentation.threshold')
            case 'zero'
                for f = tpms
                    seg.(f{1}) = (seg.(f{1}) > 0) & ~segAll;
                    segAll = segAll | seg.(f{1});
                    
                    aas_log(aap,false,sprintf('Zero thresholded image %s sums up to %d vox', ...
                        f{1}, sum(seg.(f{1})(:))))
                end
            case 'exclusive'
                % Any particular voxel has greatest chance of being...
                Ys = cellfun(@(x) seg.(x), tpms, 'UniformOutput', false);
                maxeY = max(cat(4,Ys{:}),[],4);
                
                for f = tpms
                    % We want to check where the segmentation has the
                    % greatest values
                    seg.(f{1}) = (seg.(f{1}) == maxeY) & (seg.(f{1}) > 0.01) & ~segAll;
                    segAll = segAll | seg.(f{1});
                    
                    aas_log(aap,false,sprintf('Zero thresholded image %s sums up to %d vox', ...
                        f{1}, sum(seg.(f{1})(:))))
                end
            otherwise
                thr = aas_getsetting(aap,'segmentation.threshold');
                if numel(thr) == 1, thr(1:numel(tpms)) = thr; end
                for i = 1:numel(tpms)
                    f = tpms(i);
                    maxY = max(seg.(f{1})(:));
                    if thr(i) > 1, thr(i) = prctile(seg.(f{1})(seg.(f{1})>0),thr(i),'Method','exact'); end
                    if strcmp(f{1},'air'), thr(i) = 0.5; end
                    Y = seg.(f{1}) > thr(i);
                    
                    if ~any(Y(:)) % check for bad thesholding
                        aas_log(aap, true, sprintf('ERROR: No voxels above the threshold mask (%f) [max: %f]', thr(i), maxY))
                    end
                    
                    seg.(f{1}) = Y & ~segAll;
                    segAll = segAll | seg.(f{1});
                    
                    aas_log(aap,false,sprintf('Zero thresholded image %s sums up to %d vox', ...
                        f{1}, sum(seg.(f{1})(:))))
                end
        end        
        if aas_getsetting(aap,'segmentation.scalpthreshold')
            thr = aas_getsetting(aap,'segmentation.scalpthreshold');
            if thr > 1
                try
                    thr = prctile(scalp(scalp>0),thr,'Method','exact'); 
                catch
                    thr = prctile(scalp(scalp>0),thr); 
                end
            end
            seg.softtissue = scalp > thr;
        end
        
        % - rename
        if any(strcmp(strsplit(settings.tissue,':'),'brain'))
            % - 'brain','skull','scalp'
            seg.scalp = seg.softtissue;
            seg.skull = seg.bone & ~seg.scalp; % remove overlap;
            seg.brain = (seg.gray | seg.white | seg.cfg) & ~seg.scalp; % remove overlap;
            seg = rmfield(seg,tpms);
        else
            % - 'gray','white','csf','skull','scalp'
            seg.scalp = seg.softtissue;
            seg.skull = seg.bone & ~seg.scalp; % remove overlap
            seg.csf = seg.csf & ~seg.scalp; % remove overlap
            seg.white = seg.white & ~seg.scalp; % remove overlap
            seg.gray = seg.gray & ~seg.scalp; % remove overlap
            seg = rmfield(seg,{'bone','softtissue','air'});
        end        

        % - plot
        seg_i = ft_datatype_segmentation(seg,'segmentationstyle','indexed');
        cfg              = [];
        cfg.funparameter = 'seg';
        cfg.funcolormap  = lines(6); % distinct color per tissue
        cfg.location     = 'center';
        cfg.atlas        = seg_i;    % the segmentation can also be used as atlas
        ft_sourceplot(cfg, seg_i);
        set(gcf,'position',[0,0,720 720]);
        set(gcf,'PaperPositionMode','auto');
        print(gcf,'-noui',fullfile(aas_getsubjpath(aap,subj),['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_segmentions']),'-djpeg','-r300');
        close(gcf);

        % generate headmodel
        switch aas_getsetting(aap,'method')
            case 'singlesphere'
                cfg = [];
                cfg.method = 'singlesphere';
                cfg.tissue = settings.tissue;
                headmodel  = ft_prepare_headmodel(cfg, seg);

            case 'concentricspheres'
                cfg = [];
                cfg.method = 'concentricspheres';
                cfg.tissue = strsplit(settings.tissue,':');
                cfg.order = settings.order;
                headmodel  = ft_prepare_headmodel(cfg, seg);
            
            case 'dipoli'
                cfg = [];
                cfg.method = 'dipoli';
                cfg.tissue = strsplit(settings.tissue,':');
                headmodel  = ft_prepare_headmodel(cfg, seg);
                
            case 'simbio'
                cfg        = [];
                cfg.shift  = settings.meshshift;
                cfg.method = 'hexahedral';
                cfg.downsample = settings.downsample;
                mesh = ft_prepare_mesh(cfg,seg);
                
                mesh.pos = double(mesh.pos);
                cfg        = [];
                cfg.method = 'simbio';
                cfg.tissue = strsplit(settings.tissue,':');
                cfg.conductivity = settings.conductivity;
                headmodel  = ft_prepare_headmodel(cfg, mesh);
            otherwise
                aas_log(aap,true,sprintf('method %s is not implemented',aas_getsetting(aap,'method')));
        end
       
        %% Save output
        headmodel = ft_struct2single(headmodel);
        outfn = fullfile(aas_getsubjpath(aap,subj),'headmodel.mat');
        save(outfn,'headmodel');%,'-v7.3'); % headmodel may be >2GB
        aap = aas_desc_outputs(aap,'subject',subj,'headmodel',outfn);
        
        outfn = fullfile(aas_getsubjpath(aap,subj),'segmentation.mat');
        save(outfn,'seg');
        aap = aas_desc_outputs(aap,'subject',subj,'segmentation',outfn);

        FT.rmExternal('spm12');
        FT.unload;
    case 'checkrequirements'
        if ~aas_cache_get(aap,'fieldtrip'), aas_log(aap,true,'FieldTrip is not found'); end
end
end