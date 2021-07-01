function [aap,resp]=aamod_tdt_firstlevel_statistics(aap,task,subj)
resp='';

switch task
    case 'report'
        output = strrep(aas_getstreams(aap,'output'),'permuted_','');
        fncfg = cellstr(spm_select('FPList',spm_file(aas_getfiles_bystream(aap,'subject',subj,aas_getstreams(aap,'output',1)),'path'),'.*cfg.mat$'));
        dat = load(fncfg{1});
        
        if subj == 1 % first subject
            for o = 1:numel(output)
                if  ~isfield(aap.report,sprintf('html_TDT%02d',o))
                    aap.report.(sprintf('html_TDT%02d',o)).fname = fullfile(aap.report.condir,[aap.report.fbase sprintf('_TDT%02d.htm',o)]);
                    aap = aas_report_add(aap,'C00',...
                        sprintf('<a href="%s" target=_top>%s</a><br>',...
                        aap.report.(sprintf('html_TDT%02d',o)).fname,...
                        [strjoin(unique(dat.cfg.files.labelname,'stable'),' vs ') '; Measure: ' output{o}]));
                    aap = aas_report_add(aap,sprintf('TDT%02d',o),['HEAD=' strjoin(unique(dat.cfg.files.labelname,'stable'),' vs ') '; Measure: ' output{o}]);
                end
                if ~isempty(aap.tasklist.currenttask.extraparameters.aap.directory_conventions.analysisid_suffix)
                    aap = aas_report_add(aap,sprintf('TDT%02d',o),sprintf('<h2>Branch: %s</h2>',...
                        aap.tasklist.currenttask.extraparameters.aap.directory_conventions.analysisid_suffix(2:end)));
                end
            end
        end
         
        for o = 1:numel(output)            
            taskName = [strjoin(unique(dat.cfg.files.labelname,'stable'),' vs ') '; Measure: ' output{o}];
            
            aap = aas_report_add(aap,sprintf('TDT%02d',o),['Subject: ' basename(aas_getsubjpath(aap,subj)) '<br>']);
            
            aap = aas_report_add(aap,subj,sprintf('<h4>%02d. %s</h4>',o,taskName));
            
            imgfname = fullfile(aas_getsubjpath(aap,subj),...
                sprintf('diagnostic_aamod_tdt_firstlevel_statistics_%s_overlay_3_001.jpg',output{o}));
            if exist(imgfname,'file')
                tstat = dlmread(strrep(imgfname,'_overlay_3_001.jpg','.txt'));
                
                % add image and report to subject
                aap = aas_report_add(aap, subj,'<table><tr>');
                aap = aas_report_add(aap, subj, sprintf('T = %2.2f - %2.2f</tr><tr>', tstat(1), tstat(2)));
                aap = aas_report_addimage(aap, subj, imgfname);
                aap = aas_report_add(aap,subj,'</tr></table>');
                
                % add image and report to summary
                aap = aas_report_add(aap,sprintf('TDT%02d',o),'<table><tr>');
                aap = aas_report_add(aap,sprintf('TDT%02d',o),sprintf('T = %2.2f - %2.2f</tr><tr>', tstat(1), tstat(2)));
                aap = aas_report_addimage(aap,sprintf('TDT%02d',o), imgfname);
                aap = aas_report_add(aap,sprintf('TDT%02d',o),'</tr></table>');
            end
        end
    case 'doit'
        [~, TDT] = aas_cache_get(aap,'tdt');
        TDT.load;
        
        dat = load(aas_getfiles_bystream(aap,'subject',subj,'settings'));
        bgfname = aas_getfiles_bystream(aap,'subject',subj,'structural');
        
        %% Permutation
        org_cfg = dat.cfg; % keeping the unpermuted cfg to copy parameters below
        cfg = org_cfg; % initialize new cfg like the original
        
        % update file locations
        cfg.files.mask = {aas_getfiles_bystream(aap,'subject',subj,'mask')};
        beta_loc = spm_file(aas_getfiles_bystream(aap,'subject',subj,'firstlevel_betas'),'path'); beta_loc = beta_loc(1,:);
        cfg.files.name = spm_file(cfg.files.name,'path',beta_loc);        
        
        cfg = rmfield(cfg,'design');
        cfg.design.function = org_cfg.design.function;
        
        cfg.results = rmfield(cfg.results, 'resultsname'); % the name should be generated later
        cfg.results.dir = fullfile(aas_getsubjpath(aap,subj), 'decoding', 'perm'); % change directory
        cfg.results.write = 2; % write only mat file
        cfg.results.overwrite = 1;
        
%         combine = 0;   % see make_design_permutations how you can run all analysis in one go, might be faster but takes more memory
%         designs = make_design_permutation(cfg,aas_getsetting(aap,'permutation.iteration'),combine);
%         
%         doParallel = aas_getsetting(aap,'permutation.numberofworkers') > 1;
%         if doParallel
%             aapoolprofile = strsplit(aap.directory_conventions.poolprofile,':'); poolprofile = aapoolprofile{1};
%             if ~strcmp(aap.options.wheretoprocess,'qsub'), aas_log(aap,false,sprintf('WARNING: pool profile %s is not used via DCS/MPaS; therefore it may not work for parfor',poolprofile)); end
%             try
%                 cluster = parcluster(poolprofile);
%                 if numel(aapoolprofile) > 1, cluster.ResourceTemplate = strjoin({aapoolprofile{2} cluster.ResourceTemplate}, ' '); end
%                 global aaworker;
%                 wdir = spm_file(tempname,'basename'); wdir = fullfile(aaworker.parmpath,wdir(1:8));
%                 aas_makedir(aap,wdir);
%                 cluster.JobStorageLocation = wdir;
%                 nWorkers = min([cluster.NumWorkers aas_getsetting(aap,'permutation.numberofworkers')]);
%                 pool = gcp('nocreate');
%                 if isempty(pool) || pool.NumWorkers < nWorkers, delete(pool); pool = parpool(cluster,nWorkers); end
%             catch E
%                 aas_log(aap,false,['WARNING: ' poolprofile ' could not been initialised - ' E.message ' --> parallelisation is disabled']);
%                 doParallel = false;
%             end
%             
%             if doParallel
%                 parfor i_perm = 1:aas_getsetting(aap,'permutation.iteration')
%                     cfg_perm = cfg;
%                     cfg_perm.design = designs{i_perm};
%                     cfg_perm.results.filestart = ['perm' sprintf('%04d',i_perm)];
%                     
%                     decoding(cfg_perm); % run permutation
%                 end
%                 
%                 delete(pool);
%             end
%         end
%         if ~doParallel
%             for i_perm = 1:aas_getsetting(aap,'permutation.iteration')
%                 cfg_perm = cfg;
%                 cfg_perm.design = designs{i_perm};
%                 cfg_perm.results.filestart = ['perm' sprintf('%04d',i_perm)];
%                 
%                 decoding(cfg_perm); % run permutation
%             end
%         end
        
        % - write 4D NIfTI
        dim = [cfg.datainfo.dim aas_getsetting(aap,'permutation.iteration')];
        for o = 1:numel(cfg.results.output)
            outpermfnames = cellstr(spm_select('FPList',cfg.results.dir,['^perm[0-9]{4}_' cfg.results.output{o} '.mat']));
            outfname = fullfile(cfg.results.dir,['perm_' cfg.results.output{o} '.nii']);
            
            N      = nifti;
            N.dat  = file_array(outfname,dim,[16 0]);
            N.mat  = cfg.datainfo.mat;
            N.mat0 = cfg.datainfo.mat;
            N.descrip     = ['permutation ' cfg.results.output{o} ' 4D'];
            create(N);
            
            for n = 1:aas_getsetting(aap,'permutation.iteration')
                dat = load(outpermfnames{n}); results = dat.results;
                output = results.(cfg.results.output{o}).output;
                Y = zeros(results.datainfo.dim);
                if strcmp(cfg.analysis,'roi')
                    roiInd = cellfun(@(rname) sscanf(rname,'roi%05d'), results.roi_names);
                    output = num2cell(output);
                else
                    roiInd = 1;
                    output = {output};
                end
                for r = 1:numel(roiInd)
                    Y(results.mask_index_each{r}) = output{r};
                end                
                N.dat(:,:,:,n) = Y;
                spm_get_space([N.dat.fname ',' num2str(n)], N.mat);
            end
            
            N.dat = reshape(N.dat,dim);
            
            aap = aas_desc_outputs(aap,'subject',subj,['permuted_' cfg.results.output{o}],outfname);
        end
        
        %% Statistics
        statdir = spm_file(cfg.results.dir,'basename','stats');
        aas_makedir(aap,statdir);
        for o = 1:numel(cfg.results.output)
            result = cellstr(aas_getfiles_bystream(aap,'subject',subj,cfg.results.output{o}));
            dat = load(result{strcmp(spm_file(result,'ext'),'mat')}); res_mat = dat.results;
            reference = spm_select('FPList',spm_file(aas_getfiles_bystream(aap,'subject',subj,['permuted_' cfg.results.output{o}],'output'),'path'),...
              ['^perm[0-9]{4}_' cfg.results.output{o} '.mat$']);%             
            dat = load(aas_getfiles_bystream(aap,'subject',subj,'settings')); statcfg = dat.cfg;
            statcfg.stats.test = 'permutation';
            statcfg.stats.tail = 'right';
            statcfg.stats.output = statcfg.results.output{o};
            statcfg.stats.results.write = 1;
            statcfg.stats.results.fpath = statdir;
            decoding_statistics(statcfg,res_mat,reference);
            
            % - create thresholded image
            fnstatres = fullfile(statdir,['stats_' statcfg.stats.output '_' statcfg.stats.test '_' statcfg.stats.tail '.mat']);
            statres = load(fnstatres); statres = statres.results_out;

            outfname = spm_file(fnstatres,'ext','nii','prefix','thresh_');
            N      = nifti;
            N.dat  = file_array(outfname,statres.datainfo.dim,[16 0]);
            N.mat  = statres.datainfo.mat;
            N.mat0 = statres.datainfo.mat;
            N.descrip     = ['thresholded ' cfg.results.output{o}];
            create(N);
            
            Z = statres.(cfg.results.output{o}).output;
            Z(statres.(cfg.results.output{o}).p >= 0.05) = 0;
            Y = zeros(statres.datainfo.dim);
            if strcmp(cfg.analysis,'roi')
                roiInd = cellfun(@(rname) sscanf(rname,'roi%05d'), statres.roi_names);
                Z = num2cell(Z);
            else
                roiInd = 1;
                Z = {Z};
            end
            for r = 1:numel(roiInd)
                Y(statres.mask_index_each{r}) = Z{r};
            end
            N.dat(:,:,:) = Y;
            spm_get_space(N.dat.fname, N.mat);
            
            N.dat = reshape(N.dat,statres.datainfo.dim);
            
            % - overlay
            % -- edges of activation
            slims = ones(4,2);
            sAct = arrayfun(@(x) any(Y(x,:,:),'all'), 1:size(Y,1));
            if numel(find(sAct))<2, slims(1,:) = [1 size(Y,1)];
            else, slims(1,:) = [find(sAct,1,'first') find(sAct,1,'last')]; end
            sAct = arrayfun(@(y) any(Y(:,y,:),'all'), 1:size(Y,2));
            if numel(find(sAct))<2, slims(2,:) = [1 size(Y,2)];
            else, slims(2,:) = [find(sAct,1,'first') find(sAct,1,'last')]; end
            sAct = arrayfun(@(z) any(Y(:,:,z),'all'), 1:size(Y,3));
            if numel(find(sAct))<2, slims(3,:) = [1 size(Y,3)];
            else, slims(3,:) = [find(sAct,1,'first') find(sAct,1,'last')]; end
            % -- convert to mm
            slims = sort(N.mat*slims,2);
            % -- extend if too narrow (min. 50mm)
            slims = slims + (repmat([-25 25],4,1).*repmat(diff(slims,[],2)<50,1,2));
            
            % - draw
            axis = {'sagittal','coronal','axial'};
            for a = 1:3
                if any(cell2mat(Z)~=0), stat_fname = {N.dat.fname}; else, stat_fname = {}; end
                [fig, v] = map_overlay(bgfname,stat_fname,axis{a},slims(a,1):aas_getsetting(aap,'overlay.nth_slice'):slims(a,2));
                fnsl{a} = fullfile(aas_getsubjpath(aap,subj), sprintf('diagnostic_aamod_tdt_firstlevel_statistics_%s_overlay_%d.jpg',cfg.results.output{o},a));
                
                if (~any(cell2mat(Z)~=0))
                    annotation('textbox',[0 0.475 0.5 0.5],'String','No voxels survive threshold','FitBoxToText','on','fontweight','bold','color','y','fontsize',18,'backgroundcolor','k');
                end
                
                spm_print(fnsl{a},fig,'jpg')
                close(fig);
            end
            
            dlmwrite(fullfile(aas_getsubjpath(aap,subj), sprintf('diagnostic_aamod_tdt_firstlevel_statistics_%s.txt',cfg.results.output{o})),[min(v(v~=0)), max(v)]);
        end
        
        %% Cleanup
        TDT.unload;
        
    case 'checkrequirements'
        if ~aas_cache_get(aap,'tdt'), aas_log(aap,true,'TDT is not found'); end
        
        % automatic connection to the nearest aamod_tdt_decode
        src = aas_getstreams(aas_setcurrenttask(aap,aas_getsourcestage(aap,'aamod_tdt_decode','settings')),'output'); src(1:2) = []; % settings, mask
        inp = aas_getstreams(aap,'input'); inp(1:4) = []; % settings, mask, firstlevel_betas, structural
        for s = 1:numel(src)
            if strcmp(inp{s},src{s}), continue; end
            if s == 1
                aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'input',src{s},'input');
                aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'permuted_output',['permuted_' src{s}],'output');
            else
                aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'append',src{s},'input');
                aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'append',['permuted_' src{s}],'output');
            end
            aas_log(aap,false,['INFO: ' aap.tasklist.currenttask.name ' input/output streams: ''' src{s} '''']);
        end
end
end