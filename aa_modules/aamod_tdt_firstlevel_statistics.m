function [aap,resp]=aamod_tdt_firstlevel_statistics(aap,task,subj)
resp='';

switch task
    case 'report'
        output = strrep(aas_getstreams(aap,'output'),'permuted_','');
        fncfg = cellstr(spm_select('FPList',spm_file(aas_getfiles_bystream(aap,'subject',subj,aas_getstreams(aap,'output',1)),'path'),'.*cfg.mat$'));
        dat = load(fncfg{1});
        
        if subj == 1 % first subject -> setup page
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
        %% Init
        [~, TDT] = aas_cache_get(aap,'tdt');
        TDT.load;
        
        if aas_stream_has_contents(aap,'subject',subj,'structural')
            bgfname = aas_getfiles_bystream(aap,'subject',subj,'structural');   
        else
            bgfname = '';
        end        
        load(aas_getfiles_bystream(aap,'subject',subj,'settings'),'cfg');
        cfgs = cfg;
        dopairwise = aas_stream_has_contents(aap,'subject',subj,'settings_pairwise');
        npw = 0;
        if dopairwise            
            for fnpw = cellstr(aas_getfiles_bystream(aap,'subject',subj,'settings_pairwise'))'
                load(fnpw{1},'cfg');
                npw = npw + 1;
                cfgs(1+npw) = cfg;
            end
        end
        doonevsall = aas_stream_has_contents(aap,'subject',subj,'settings_onevsall');
        nova = 0;
        if doonevsall            
            for fnova = cellstr(aas_getfiles_bystream(aap,'subject',subj,'settings_onevsall'))'
                load(fnova{1},'cfg');
                nova = nova + 1;
                cfgs(1+npw+nova) = cfg;
            end
        end
        
        if aas_getsetting(aap,'permutation.iteration') >= 1e4
            aas_log(aap,false,'WARNING: permutation with more than 9999 iterations is not supported -> permutation.iteration = 9999');
            aap.tasklist.currenttask.settings.permutation.iteration = 9999;
        end
        
        srcmodulename = strrep(aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name,'firstlevel_statistics','decode'); % CAVE: assumption on module name
        srcpath = aas_getsubjpath(aas_setcurrenttask(aap,aas_getsourcestage(aap,srcmodulename,'settings')),subj);
        
        %% Generate prescription (cfg, design, doIter)
        toExclude = [];
        doIter = cell(1,numel(cfgs));
        for indCfg = 1:numel(cfgs)
            org_cfg = cfgs(indCfg); % keeping the unpermuted cfg to copy parameters below
            
            % update file locations
            cfgs(indCfg).files.mask = {aas_getfiles_bystream(aap,'subject',subj,'mask')};
            cfgs(indCfg).files.name = strrep(cfgs(indCfg).files.name,srcpath,aas_getsubjpath(aap,subj));
            
            cfgs(indCfg).design = [];
            cfgs(indCfg).design.function = org_cfg.design.function;
            
            cfgs(indCfg).results = rmfield(cfgs(indCfg).results, 'resultsname'); % the name should be generated later
            cfgs(indCfg).results.dir = fullfile(aas_getsubjpath(aap,subj), 'decoding', 'perm',cfgs(indCfg).results.filestart); % change directory
            cfgs(indCfg).results.write = 2; % write only mat file
            cfgs(indCfg).results.overwrite = 1;
            
            aas_makedir(aap,cfgs(indCfg).results.dir);
            
            fnDesign = fullfile(cfgs(indCfg).results.dir,[cfgs(indCfg).results.filestart '_designs.mat']);
            
            if exist(fnDesign,'file') && isempty(spm_select('List',cfgs(indCfg).results.dir,['^' cfgs(indCfg).results.filestart '_perm[0-9]{4}_cfg.mat']))
                delete(fnDesign);                
            end
            
            if exist(fnDesign,'file')
                clear designs; load(fnDesign,'designs');
                doIter{indCfg} = 1:numel(designs);
                doneIter = cellfun(@(f) sscanf(f,[cfgs(indCfg).results.filestart '_perm%04d']), cellstr(spm_select('List',cfgs(indCfg).results.dir,'.*_cfg.mat')));
                doIter{indCfg} = setdiff(doIter{indCfg}, doneIter);
            else
                combine = 0;   % see make_design_permutations how you can run all analysis in one go, might be faster but takes more memory
                designs = make_design_permutation(cfgs(indCfg),aas_getsetting(aap,'permutation.iteration'),combine);
                save(fnDesign,'designs')
                doIter{indCfg} = 1:numel(designs);
            end
            if numel(designs) < aas_getsetting(aap,'permutation.iteration')
                aas_log(aap,false,sprintf('WARNING: Number of possible iterations (%d) is lower than requested (%d) -> no permutations will be performed.',numel(designs),aas_getsetting(aap,'permutation.iteration')));
                toExclude(end+1) = indCfg;
                rmdir(cfgs(indCfg).results.dir,'s');
            end
        end
        cfgs(toExclude) = []; doIter(toExclude) = [];
        
        %% Run permutation
        doParallel = aas_getsetting(aap,'permutation.numberofworkers') > 1;
        if doParallel
            aapoolprofile = strsplit(aap.directory_conventions.poolprofile,':'); poolprofile = aapoolprofile{1};
            if ~strcmp(aap.options.wheretoprocess,'batch'), aas_log(aap,false,sprintf('WARNING: pool profile %s is not used via DCS/MPaS; therefore it may not work',poolprofile)); end
            try
                global aaparallel;
                aaparallel = aap.options.aaparallel;
                global aaworker;
                q = aaq_batch(aap);
                cluster = q.pool;
                
                wdir = spm_file(tempname,'basename'); wdir = fullfile(aaworker.parmpath,wdir(1:8));
                aas_log(aap,false,['INFO: parallel job storage is ' wdir]);
                aas_makedir(aap,wdir);
                cluster.JobStorageLocation = wdir;
                if cluster.NumWorkers > 1
                    nWorkers = min([cluster.NumWorkers aas_getsetting(aap,'permutation.numberofworkers')]);
                else
                    nWorkers = aas_getsetting(aap,'permutation.numberofworkers');
                end
                aas_log(aap,false,sprintf('INFO: %s initialised with %d workers.',poolprofile,nWorkers));
            catch E
                aas_log(aap,false,['WARNING: ' poolprofile ' could not been initialised - ' E.message ' --> parallelisation is disabled']);
                doParallel = false;
            end
        end
        if doParallel
            maxIterPerJob = ceil(sum(cellfun(@numel, doIter))/nWorkers); 
            cfg_perm = cfgs(1); cfg_perm(1) = []; % create empty array of cfgs
            for indCfg = 1:numel(cfgs)
                clear designs; load(fullfile(cfgs(indCfg).results.dir,[cfgs(indCfg).results.filestart '_designs.mat']),'designs');
                for i_perm = doIter{indCfg}                    
                    cfg_perm(end+1) = cfgs(indCfg);
                    cfg_perm(end).design = designs{i_perm};
                    if indCfg > 1+npw, cfg_perm(end).design.unbalanced_data = 'ok'; end
                    cfg_perm(end).results.filestart = [cfgs(indCfg).results.filestart '_perm' sprintf('%04d',i_perm)];

                    if (numel(cfg_perm) == maxIterPerJob) ||... % batch is ready
                            (indCfg == numel(cfgs) && i_perm == max(doIter{indCfg})) % last batch
                        batch(cluster,@batch_decoding,0,{cfg_perm},'AutoAttachFiles',false); % run permutations
                        cfg_perm(1:end) = []; % reset array of cfgs
                    end
                end
            end
            while ~all(arrayfun(@(j) strcmp(j.State,'finished'), cluster.Jobs)), pause(1); end
        else
            for indCfg = 1:numel(cfgs)
                clear designs; load(fullfile(cfgs(indCfg).results.dir,[cfgs(indCfg).results.filestart '_designs.mat']),'designs');
                for i_perm = doIter{indCfg}                    
                    cfg_perm = cfgs(indCfg);
                    cfg_perm.design = designs{i_perm};
                    if indCfg > 1+npw, cfg_perm(end).design.unbalanced_data = 'ok'; end
                    cfg_perm.results.filestart = [cfgs(indCfg).results.filestart '_perm' sprintf('%04d',i_perm)];
                    
                    decoding(cfg_perm); % run permutation
                end
            end
        end
        
        %% Write 4D NIfTI
        fnpw = {}; fnova = {};
        for indCfg = 1:numel(cfgs)
            clear designs; load(fullfile(cfgs(indCfg).results.dir,[cfgs(indCfg).results.filestart '_designs.mat']),'designs');
            dim = [cfgs(indCfg).datainfo.dim numel(designs)];
            for o = 1:numel(cfgs(indCfg).results.output)
                aas_log(aap,false,['INFO: Saving - ' spm_file(cfgs(indCfg).results.dir,'basename') ' ' cfgs(indCfg).results.output{o}]);
                outpermfnames = cellstr(spm_select('FPList',cfgs(indCfg).results.dir,['^' cfgs(indCfg).results.filestart '_perm[0-9]{4}_' cfgs(indCfg).results.output{o} '.mat']));
                outfname = fullfile(cfgs(indCfg).results.dir,[cfgs(indCfg).results.filestart '_perm_' cfgs(indCfg).results.output{o} '.nii']);
                
                N      = nifti;
                N.dat  = file_array(outfname,dim,[16 0]);
                N.mat  = cfgs(indCfg).datainfo.mat;
                N.mat0 = cfgs(indCfg).datainfo.mat;
                N.descrip     = ['permutation ' cfgs(indCfg).results.output{o} ' 4D'];
                create(N);
                
                for n = 1:numel(designs)
                    dat = load(outpermfnames{n}); results = dat.results;
                    output = results.(cfgs(indCfg).results.output{o}).output;
                    Y = zeros(results.datainfo.dim);
                    if strcmpi(cfgs(indCfg).analysis,'roi') && (numel(results.roi_names) > 1)
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
                
                if indCfg == 1, aap = aas_desc_outputs(aap,'subject',subj,['permuted_' cfgs(indCfg).results.output{o}],outfname); 
                elseif indCfg <= 1+npw, fnpw(end+1) = cellstr(outfname);
                elseif indCfg <= 1+npw+nova, fnova(end+1) = cellstr(outfname);
                end
            end
            clear N
        end
        if ~isempty(fnpw), aap = aas_desc_outputs(aap,'subject',subj,['permuted_' cfgs(indCfg).results.output{o} '_pairwise'],fnpw); end
        if ~isempty(fnova), aap = aas_desc_outputs(aap,'subject',subj,['permuted_' cfgs(indCfg).results.output{o} '_onevsall'],fnova); end
        
        %% Statistics
        if aas_getsetting(aap,'dostatistics')
            cfg = cfgs(1); % for main analysis only (TODO: add for the pairwise, if they exist)
            statdir = spm_file(cfg.results.dir,'basename','stats');
            aas_makedir(aap,statdir);
            for o = 1:numel(cfg.results.output)
                aas_log(aap,false,['INFO: Inferencing - ' cfg.results.output{o}]);
                result = cellstr(aas_getfiles_bystream(aap,'subject',subj,cfg.results.output{o}));
                dat = load(result{strcmp(spm_file(result,'ext'),'mat') & strcmp(spm_file(result,'basename'),spm_file(cfg.results.output{o},'prefix','res_'))}); res_mat = dat.results; % TODO: potential issue with connectivity, where both results are MAT files
                reference = spm_select('FPList',spm_file(aas_getfiles_bystream(aap,'subject',subj,['permuted_' cfg.results.output{o}],'output'),'path'),...
                    ['^res_perm[0-9]{4}_' cfg.results.output{o} '.mat$']);%
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
                V.dt = [spm_type('float32') 0];
                V.mat = statres.datainfo.mat;
                V.private.mat = statres.datainfo.mat;
                
                Z = statres.(cfg.results.output{o}).z;
                thrZ = Z; thrZ(statres.(cfg.results.output{o}).p >= 0.05) = 0;
                if strcmp(cfg.analysis,'roi') && (numel(results.roi_names) > 1)
                    roiInd = cellfun(@(rname) sscanf(rname,'roi%05d'), statres.roi_names);
                    Z = num2cell(Z);
                    thrZ = num2cell(thrZ);
                else
                    roiInd = 1;
                    Z = {Z};
                    thrZ = {thrZ};
                end
                
                Y = zeros(statres.datainfo.dim);
                for r = 1:numel(roiInd)
                    Y(statres.mask_index_each{r}) = Z{r};
                end
                outfname = spm_file(fnstatres,'ext','nii');
                nifti_write(outfname,Y,cfg.results.output{o},V);
                aap = aas_desc_outputs(aap,'subject',subj,['zstat_' cfg.results.output{o}],outfname);

                Y = zeros(statres.datainfo.dim);
                for r = 1:numel(roiInd)
                    Y(statres.mask_index_each{r}) = thrZ{r};
                end
                outfname = spm_file(fnstatres,'ext','nii','prefix','thresh_');
                nifti_write(outfname,Y,['thresholded ' cfg.results.output{o}],V);
                aap = aas_desc_outputs(aap,'subject',subj,['thresholded_' cfg.results.output{o}],outfname);
                
                % - overlay
                if ~isempty(bgfname)
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
                    slims = sort(V.mat*slims,2);
                    % -- extend if too narrow (min. 50mm)
                    slims = slims + (repmat([-25 25],4,1).*repmat(diff(slims,[],2)<50,1,2));
                    
                    % - draw
                    axis = {'sagittal','coronal','axial'};
                    for a = 1:3
                        if any(cell2mat(Z)~=0), stat_fname = {outfname}; else, stat_fname = {}; end
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
            end
        end
        
        %% Cleanup
        TDT.unload;
        
    case 'checkrequirements'
        if subj == 1 % only once
            if ~aas_cache_get(aap,'tdt'), aas_log(aap,true,'TDT is not found'); end
            
            % automatic connection to the nearest aamod_tdt_decode
            srcmodulename = strrep(aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name,'firstlevel_statistics','decode'); % CAVE: assumption on module name
            src = aas_getstreams(aas_setcurrenttask(aap,aas_getsourcestage(aap,srcmodulename,'settings')),'output'); src = setdiff(src,{'settings' 'samples' 'mask'});
            inp = aas_getstreams(aap,'input');
            if any(strcmp(src,'settings_pairwise'))
                if ~any(strcmp(inp,'settings_pairwise')), aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'append','settings_pairwise','input'); end
                src = setdiff(src,{'settings_pairwise'});
            end
            if any(strcmp(src,'settings_onevsall'))
                if ~any(strcmp(inp,'settings_onevsall')), aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'append','settings_onevsall','input'); end
                src = setdiff(src,{'settings_onevsall'});
            end
            inp(contains(inp,{'settings' 'settings_pairwise' 'settings_onevsall' 'samples' 'mask' 'firstlevel_betas' 'roidata_firstlevel_betas' 'valid_roi_firstlevel_betas' 'connectivity' 'structural'})) = [];
            for s = 1:numel(src)
                if s <= numel(inp) && strcmp(inp{s},src{s}), continue; end
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
end

function batch_decoding(cfgs)
for c = reshape(cfgs,1,[])
    decoding(c)
end
end