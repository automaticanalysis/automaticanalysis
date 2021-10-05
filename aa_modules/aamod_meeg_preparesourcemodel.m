function [aap, resp] = aamod_meeg_preparesourcemodel(aap,task,subj)

resp='';

switch task
    case 'report'
        
    case 'doit'
        [~, FT] = aas_cache_get(aap,'fieldtrip');
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
                tmpltfile = fullfile(FT.toolPath, ['template/sourcemodel/standard_sourcemodel3d' num2str(aas_getsetting(aap,'options.grid.resolution')) 'mm.mat']);
                dat = load(tmpltfile);
                cfg           = [];
                cfg.warpmni   = 'yes';
                cfg.spmmethod = 'new';
                cfg.template  = dat.sourcemodel;
                cfg.nonlinear = 'yes';
                cfg.mri       = mri;
                cfg.unit      ='mm';
                sourcemodel   = ft_prepare_sourcemodel(cfg);
                sourcemodel.cfg = [];
                sourcemodel.cfg.template = tmpltfile;
                sourcemodel.cfg.urinside = sourcemodel.inside;
                sourcemodel.pos(~sourcemodel.inside,:) = [];
                sourcemodel.inside(~sourcemodel.inside,:) = [];
                
            case 'corticalsheet'
                aas_runFScommand(aap,sprintf('export PATH=$PATH:%s/bin_rh_linux64; %s/bin/ft_postfreesurferscript.sh %s %s %s',WB.toolPath, FT.toolPath,...
                    aas_getstudypath(aap),aas_getsubjname(aap,subj),...
                    fullfile(WB.templateDir,'standard_mesh_atlases')));
                resPerHemi = sscanf(aas_getsetting(aap,'options.corticalsheet.resolution'),'%dk')/2; downsample = 1;
                if resPerHemi < 4 % create 4k meshes and downsample
                    downsample = 4/resPerHemi;
                    resPerHemi = 4;
                end
                sheetfn = fullfile(aas_getsubjpath(aap,subj),'workbench',sprintf('%s.L.midthickness.%dk_fs_LR.surf.gii',aas_getsubjname(aap,subj),resPerHemi));
                if ~exist(sheetfn,'file'), aas_log(aap,true,'ft_postfreesurferscript failed'); end
                sourcemodel = ft_read_headshape({sheetfn, strrep(sheetfn, '.L.', '.R.')});
                sourcemodel.inside = sourcemodel.atlasroi>0;
                sourcemodel = rmfield(sourcemodel, 'atlasroi');
                if downsample > 1 % run downsample if needed
                    indpos = true(1,size(sourcemodel.pos,1)); 
                    indpos(downsample:downsample:end) = false;
                    
                    origpos = sourcemodel.pos;
                    for f = intersect(fieldnames(sourcemodel),{'pos' 'inside' 'sulc' 'curv' 'thickness' 'brainstructure'})'
                        sourcemodel.(f{1})(indpos,:) = [];
                    end
                    
                    % process tri
                    % - resample nodes to nearest
                    tri = arrayfun(@(x) dsearchn(sourcemodel.pos, origpos(x,:)),sourcemodel.tri(:));
                    tri = reshape(tri,[],3);
                    % - remove duplicates
                    [~,ia] = unique(sort(tri,2),'stable','rows');
                    tri = tri(ia,:);
                    indtri = arrayfun(@(t) numel(unique(tri(t,:))) < 3, 1:size(tri,1));
                    tri(indtri,:) = [];                    
                    sourcemodel.tri = tri;
                end
                
                % write surfaces
                fnames = {};
                cfg = [];
                cfg.parameter = 'param';
                cfg.precision = 'single';
                for surftype = {'inflated' 'midthickness' 'pial','white'}
                    sheetfn = fullfile(aas_getsubjpath(aap,subj),'workbench',sprintf('%s.L.%s.164k_fs_LR.surf.gii',aas_getsubjname(aap,subj),surftype{1}));
                    surf = ft_read_headshape({sheetfn, strrep(sheetfn, '.L.', '.R.')});
                    surf.(cfg.parameter) = ones(size(surf.pos,1),1);
                    cfg.filename = fullfile(aas_getsubjpath(aap,subj),['surf_' surftype{1}]);
                    ft_sourcewrite(cfg,surf);
                    fnames = vertcat(fnames,{[cfg.filename '.gii']});
                end
                aap = aas_desc_outputs(aap,'subject',subj,'sourcesurface',fnames);
                
                % create atlas
                fnames = cellstr(aas_getfiles_bystream(aap,'subject',subj,'freesurfer'));
                res = regexp(fnames,['.*' aas_getsetting(aap,'options.corticalsheet.annotation') '\.annot$'],'match');
                annot = sort(vertcat(res{:}));
                res = regexp(fnames,'.*h\.pial$','match');
                mesh = sort(vertcat(res{:}));
                
                atlaslh = ft_read_atlas({annot{1} mesh{1}},'format','freesurfer_aparc');
                atlasrh = ft_read_atlas({annot{2} mesh{2}},'format','freesurfer_aparc');
                
                sourceatlas = rmfield(atlaslh,'rgba');
                sourceatlas.pos = vertcat(atlaslh.pos, atlasrh.pos);
                sourceatlas.tri = vertcat(atlaslh.tri, atlasrh.tri+size(atlaslh.pos,1));
                sourceatlas.aparclabel = vertcat(spm_file(atlaslh.aparclabel,'prefix','lh_'),spm_file(atlaslh.aparclabel,'prefix','rh_'));
                sourceatlas.aparc = vertcat(atlaslh.aparc, atlasrh.aparc+numel(atlaslh.aparclabel));
                
                % diag
                fname = fullfile(aas_getsubjpath(aap,subj),['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_atlas.jpg']);
                if ~exist(fname,'file')
                    cfg = [];
                    cfg.figure = figure; hold on;
                    cfg.method = 'surface';
                    cfg.funparameter = 'aparc';
                    cfg.funcolormap = distinguishable_colors(numel(sourceatlas.aparclabel));
                    ft_sourceplot(cfg,sourceatlas);
                    view(135,45);
                    set(cfg.figure,'position',[0,0,720 720]);
                    set(cfg.figure,'PaperPositionMode','auto');
                    print(cfg.figure,'-noui',fname,'-djpeg','-r300');
                    close(cfg.figure);
                end
                
                fnAtlas = fullfile(aas_getsubjpath(aap,subj),'sourceatlas.mat');
                save(fnAtlas,'sourceatlas');
                aap = aas_desc_outputs(aap,'subject',subj,'sourceatlas',fnAtlas);
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
        if subj == 1
            [s, FT] = aas_cache_get(aap,'fieldtrip');
            if ~s, aas_log(aap,true,'FieldTrip is not found'); end
            if strcmp(aas_getsetting(aap,'method'),'corticalsheet')
                [s, WB] = aas_cache_get(aap,'hcpwb');
                if ~s, aas_log(aap,true,'HCP Workbench is not found -> corticalsheet is not available'); end
                if isempty(WB.templateDir) || ~exist(WB.templateDir,'dir'), aas_log(aap,true,'templates for HCP Workbench are not found -> corticalsheet is not available');
                else
                    aas_log(aap,false,sprintf('WARNING: make sure that the template directory is prepared as described in %s/bin/ft_postfreesurferscript.sh',FT.toolPath))
                end
            else
                aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'sourcesurface',[],'output');
            end
        end
end
end