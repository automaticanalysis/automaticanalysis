function [aap, resp] = aamod_meeg_preparesourcemodel(aap,task,subj)

resp='';

switch task
    case 'report'
        
    case 'doit'
        [junk, FT] = aas_cache_get(aap,'fieldtrip');
        FT.load;
        
        if strcmp(aas_getsetting(aap,'method'),'corticalsheet')
            [status, WB] = aas_cache_get(aap,'hcpwb');
            if status, WB.load; end
            if ~aas_stream_has_contents(aap,'freesurfer') || ~status || isempty(WB.templateDir)
                aas_log(aap,false,'Freesurfer data, HCP Workbench or template directory is not found -> corticalsheet is not available');
                aap.tasklist.currenttask.settings.method = 'grid';
            end
        end
        
        % load MRI
        if strcmp(aas_getsetting(aap,'method'),'corticalsheet')
            aas_log(aap,false,'INFO: individual sourcemodel based on cortical sheet');
            FSfiles = cellstr(aas_getfiles_bystream(aap,'subject',subj,'freesurfer'));
            mrifn = FSfiles{strcmp(spm_file(FSfiles,'filename'),'T1.mgz')};
            mri = ft_read_mri(mrifn,'dataformat','freesurfer_mgz');
        elseif aas_stream_has_contents(aap,'structural')
            aas_log(aap,false,'INFO: individual sourcemodel based on 3D grid');
            aap.tasklist.currenttask.settings.method = 'grid'; % overwrite setting if needed
            mrifn = aas_getfiles_bystream(aap,'subject',subj,'structural');
            mri = ft_read_mri(mrifn,'dataformat','nifti_spm');
        else % fall back to template [TODO]
            aas_log(aap,true,'no structural input found -> template must be used - NYI')
        end
        mri.coordsys = 'ras';
        cfg     = [];
        cfg.dim = mri.dim;
        mri     = ft_volumereslice(cfg,mri);
        
        % create sourcemodel
        switch aas_getsetting(aap,'method')
            case 'grid'
                dat = load(fullfile(FT.toolPath, ['template/sourcemodel/standard_sourcemodel3d' num2str(aas_getsetting(aap,'options.grid.resolution')) 'mm.mat']));
                cfg           = [];
                cfg.warpmni   = 'yes';
                cfg.spmmethod = 'new';
                cfg.template  = dat.sourcemodel;
                cfg.nonlinear = 'yes';
                cfg.mri       = mri;
                cfg.unit      ='mm';
                sourcemodel   = ft_prepare_sourcemodel(cfg);
                
            case 'corticalsheet'
                aas_runFScommand(aap,sprintf('export PATH=$PATH:%s/bin_rh_linux64; %s/bin/ft_postfreesurferscript.sh %s %s %s',WB.toolPath, FT.toolPath,...
                    aas_getstudypath(aap),aas_getsubjname(aap,subj),...
                    fullfile(WB.templateDir,'standard_mesh_atlases')));
                resPerHemi = sscanf(aas_getsetting(aap,'options.corticalsheet.resolution'),'%dk')/2;
                sheetfn = fullfile(aas_getsubjpath(aap,subj),'workbench',sprintf('%s.L.midthickness.%dk_fs_LR.surf.gii',aas_getsubjname(aap,subj),resPerHemi));
                if ~exist(sheetfn,'file'), aas_log(aap,true,'ft_postfreesurferscript failed'); end
                sourcemodel = ft_read_headshape({sheetfn, strrep(sheetfn, '.L.', '.R.')});
                sourcemodel.inside = sourcemodel.atlasroi>0;
                sourcemodel = rmfield(sourcemodel, 'atlasroi');                
        end
        
        sourcemodel = ft_struct2single(sourcemodel);
        outfn = fullfile(aas_getsubjpath(aap,subj),'sourcemodel.mat');
        save(outfn,'sourcemodel');
        aap = aas_desc_outputs(aap,'subject',subj,'sourcemodel',outfn);
        
        mri = ft_struct2single(mri);
        outfn = fullfile(aas_getsubjpath(aap,subj),'structural.mat');
        save(outfn,'mri');
        aap = aas_desc_outputs(aap,'subject',subj,'structural',outfn);

        FT.unload;
    case 'checkrequirements'
        [s, FT] = aas_cache_get(aap,'fieldtrip');
        if ~s, aas_log(aap,true,'FieldTrip is not found'); end
        [s, WB] = aas_cache_get(aap,'hcpwb');
        if ~s, aas_log(aap,false,'HCP Workbench is not found -> corticalsheet is not available'); end
        if isempty(WB.templateDir) || ~exist(WB.templateDir,'dir'), aas_log(aap,false,'templates for HCP Workbench are not found -> corticalsheet is not available'); 
        else
            aas_log(aap,false,sprintf('WARNING: if you want to use cotricalsheet, make sure that the template directory is prepared as described in %s/bin/ft_postfreesurferscript.sh',FT.toolPath))
        end
end
end