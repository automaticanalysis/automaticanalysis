function [aap,resp]=aamod_tdt_decode(aap,task,subj)
resp='';

switch task
    case 'report'
        %         localpath = aas_getpath_bydomain(aap,aap.tasklist.currenttask.domain,[subj,sess]);
        %
        %         fdiag = dir(fullfile(localpath,'diagnostic_*.jpg'));
        %         if isempty(fdiag)
        %             streams=aas_getstreams(aap,'output');
        %             for streamind=1:length(streams)
        %                 % obtain output
        %                 outputfnames = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,[subj sess],streams{streamind},'output');
        %
        %                 % perform diagnostics
        %                 do_diag(outputfnames);
        %             end
        %             fdiag = dir(fullfile(localpath,'diagnostic_*.jpg'));
        %         end
        %
        %         for d = 1:numel(fdiag)
        %             aap = aas_report_add(aap,subj,'<table><tr><td>');
        %             imgpath = fullfile(localpath,fdiag(d).name);
        %             aap=aas_report_addimage(aap,subj,imgpath);
        %             aap = aas_report_add(aap,subj,'</td></tr></table>');
        %         end
    case 'doit'
        [~, TDT] = aas_cache_get(aap,'tdt');
        TDT.load;
        
        flagsReslice = aap.spm.defaults.realign.write;
        flagsReslice.which = [1 0];
        flagsReslice.interp = 0;
        
        cfg = TDT.defaults;
        cfg.analysis = aas_getsetting(aap,'method');
        
        % Disable plots
        cfg.plot_selected_voxels = 0;
        cfg.plot_design = 0;
        
        % Set and configure feature selection
        settings = aas_getsetting(aap,'featureselection');
        if isstruct(settings) && ~strcmp(settings.method,'none')
            cfg.feature_selection.n_vox = settings.numberofvoxels;
            cfg.feature_selection.optimization_criterion = settings.criterion;
            switch settings.estimation
                case 'train'
                    cfg.feature_selection.estimation = 'across';
                case 'traintest'
                    cfg.feature_selection.estimation = 'all';                    
            end
            cfg.feature_selection.scale.method = settings.scaling;
            if ~strcmp(cfg.feature_selection.scale.method,'none'), cfg.feature_selection.scale.estimation = 'all'; end
            switch settings.method
                case 'maxresponse'
                    cfg.feature_selection.method = 'filter';
                    cfg.feature_selection.filter = 'F';
                case 'maxpooledredsponse'
                    cfg.feature_selection.method = 'filter';
                    cfg.feature_selection.filter = 'F0';
                case 'nonparammaxresponse'
                    cfg.feature_selection.method = 'filter';
                    cfg.feature_selection.filter = 'U';
                case 'trainedweights'
                    cfg.feature_selection.method = 'filter';
                    cfg.feature_selection.filter = 'W';
                case 'RFE'
                    cfg.feature_selection.method = 'embedded';
                    cfg.feature_selection.embedded = 'RFE';
                    cfg.feature_selection.direction = settings.options.RFE.direction;
                    cfg.feature_selection.nested_n_vox = settings.options.RFE.nestednumberofvoxels;
                otherwise % external
                    cfg.feature_selection.method = 'filter';
                    cfg.feature_selection.filter = 'external';
                    if contains(settings.method,':')
                        cfg.feature_selection.external_fname = strsplit(settings.method,':');
                    else
                        cfg.feature_selection.external_fname = cellstr(settings.method);
                    end
                    if ~all(cellfun(@(f) exist(f,'file'), cfg.feature_selection.external_fname))
                        aas_log(aap,true,'One or more external image for feature selection not found.');
                    end
            end
        end
        
        % Set and configure classifier
        % cfg.decoding.train.classification.model_parameters = '-s 1 -c 1 -q';
        cfg.decoding.software = aas_getsetting(aap,'decoding.software');        
        settings = aas_getsetting(aap,['decoding.options.' cfg.decoding.software]);
        switch cfg.decoding.software
            case 'liblinear'
                cfg.decoding.method = 'classification';
                                
                if strcmp(settings.solver,'CS-SVC'), solver = 4;
                else
                    switch settings.solver
                        case 'LR'
                            solver = [0 6 7];
                        case 'SVC'
                            solver = [1 2 3 5];
                        case 'SVR'
                            solver = [11 12 13];
                    end
                    if settings.dual
                        solver = intersect(solver,[1 3 5 6 7 12 13]);
                    else
                        solver = intersect(solver,[0 2 5 6 11]);
                    end
                    switch settings.regularisation
                        case 'L1'
                            solver = intersect(solver,[5 6]);
                        case 'L2'
                            solver = intersect(solver,[0 1 2 3 7 11 12 13]);
                    end
                    switch settings.loss
                        case 'L1'
                            solver = intersect(solver,[3 13]);
                        case 'L2'
                            solver = intersect(solver,[1 2 5 11 12]);
                        case 'none'
                            solver = intersect(solver,[0 6 7]);
                    end
                end
                cfg.decoding.train.classification.model_parameters = sprintf('-s %d -c %d -q',solver,settings.cost);
            case 'libsvm'
                schema = aap.schema.tasksettings.aamod_tdt_decode(aap.tasklist.currenttask.index).decoding.options.libsvm;
                
                if strcmp(settings.kernel,'precomputed')
                    cfg.decoding.method = 'classification_kernel';
                else
                    cfg.decoding.method = 'classification';
                end
                
                if ischar(settings.kernelparameters.gamma) % 'auto'
                    settings.kernelparameters.gamma = 1/numel(aas_getsetting(aap,'itemList'));
                end
                
                cfg.decoding.train.(cfg.decoding.method).model_parameters = sprintf('-s %d -t %d -d %d -g %1.3f -r %1.3f -c %d -b %d -q',...
                    find(strcmp(strsplit(schema.svm.ATTRIBUTE.options,'|'),settings.svm))-1, ...
                    find(strcmp(strsplit(schema.kernel.ATTRIBUTE.options,'|'),settings.kernel))-1, ...
                    settings.kernelparameters.degree, ...
                    settings.kernelparameters.gamma ,...
                    settings.kernelparameters.coef0, ...
                    settings.cost, ...
                    contains(settings.svm,'SRV') ...
                    );
        end
        
        %% Prepare
        inps = aas_getstreams(aap,'input');
        inps = inps(logical(cellfun(@(x) exist(aas_getinputstreamfilename(aap,'subject',subj,x),'file'), inps)));
        
        if aas_stream_has_contents(aap,'subject',subj,'firstlevel_spm')
            fnSPM = aas_getfiles_bystream(aap,'subject',subj,'firstlevel_spm');
            load(fnSPM,'SPM');
            dirSPM = fullfile(aas_getsubjpath(aap,subj),aap.directory_conventions.stats_singlesubj);
            orderBF = SPM.xBF.order;
            regressor_names = design_from_spm(dirSPM);
        end
        
        % check mask
        % ensure overlap with data
        if any(contains(inps,{'mask','rois'}))
            fnMask = cellfun(@(x) aas_getfiles_bystream(aap,'subject',subj,x), inps(contains(inps,{'mask','rois'})),'UniformOutput',false);
            Ydata = [];
            if aas_stream_has_contents(aap,'subject',subj,'firstlevel_betas')
                fnBeta = cellstr(aas_getfiles_bystream(aap,'subject',subj,'firstlevel_betas'));
                Ydata = spm_read_vols(spm_vol(fnBeta{1}));
            end
            if (numel(fnMask) > 1) && ~strcmp(cfg.analysis,'roi')
                brain_mask = spm_imcalc(spm_vol(char(fnMask)),fullfile(TASKROOT,'brain_mask.nii'),'min(X)',{1});
                fnMask = {brain_mask.fname};
            end
            for f = 1:numel(fnMask)
                V = spm_vol(fnMask{f});
                if ~isempty(Ydata) && ~isequal(V.dim, size(Ydata)) % reslice mask if needed
                    aas_log(aap,false,sprintf('There is a geometry mismatch between the mask (%s) and the data (%s) -> reslicing mask (assuming alignment)\n.',fnMask{f},fnBeta{1}));
                    spm_reslice({fnBeta{1} fnMask{f}},flagsReslice);
                    fnMask{f} = spm_file(fnMask{f},'prefix',flagsReslice.prefix);
                    V = spm_vol(fnMask{f});
                end
                
                Y = spm_read_vols(V);
                mask_num = unique(Y(:));
                mask_num(isnan(mask_num)|mask_num==0) = [];
                if numel(mask_num) > 1 && ~strcmp(cfg.analysis,'roi')
                    fnMask{f} = spm_file(V.fname,'prefix','b_');
                    Y = Y .* (Y > 0.5);
                end
                fnMask{f} = spm_file(fnMask{f},'prefix','ok_');
                if ~isempty(Ydata), Y(isnan(Ydata)) = 0; end
                V.fname = fnMask{f};
                V.pinfo = [1 0 0]';
                spm_write_vol(V,Y);
            end
            
            cfg.files.mask = fnMask;
        end
        
        labelnames = aas_getsetting(aap,'itemList');
        
        if aas_stream_has_contents(aap,'subject',subj,'firstlevel_betas') % activity
            % check labelnames for SPM design
            for l = 1:numel(labelnames)
                if numel(labelnames{l}) > 1, labelnames{l} = ['regexp:^[' sprintf('(%s)',labelnames{l}{:}) ']*$']; end
                if orderBF > 1 && isempty(regexp(labelnames{l},'(?<=bin)[ \[\(]{1,3}[1-9]','once'))
                    if ~isempty(strfind(labelnames{l},'regexp')), aas_log(aap,true,'Regular expression and multiple event names do not support automatic detection of multi-order BF'); end
                    labelnames{l} = ['regexp:^' strrep(labelnames{l},'*','.*') 'bin [' sprintf('(%d)',1:9) ']']; % add all orders as defaults
                end
                if iscell(labelnames{l}), labelnames{l} = labelnames{l}{1}; end
            end
            
            switch cfg.analysis
                case 'searchlight'
                    cfg.searchlight = aas_getsetting(aap,'searchlight');
                case 'network'
                    cfg.analysis = 'roi';
                    netspec = aas_getsetting(aap,'network.networkrois');
                    if ~isempty(netspec) % single-volume multi-label roi expected
                        aas_log(aap,false,'INFO: recoding rois for networks')
                        V = spm_vol(cfg.files.mask{1});
                        Y = spm_read_vols(V);
                        M = zeros(size(Y));
                        for net = 1:numel(netspec)
                            m = false(size(Y));
                            for r = reshape(netspec{net},1,[])
                                m = m | (Y == r);
                            end
                            M = M + net*m;
                        end
                    end
                    spm_write_vol(V,M);
            end            
        else
            cfg.software = 'matlab';
            cfg.analysis = 'roi';

            if aas_stream_has_contents(aap,'subject',subj,'roidata_firstlevel_betas') % roi
                labelnamesROI = labelnames;
                if isstruct(labelnamesROI), labelnamesROI = labelnamesROI.roi; end
                % check labelnames for SPM design
                for l = 1:numel(labelnamesROI)
                    if numel(labelnamesROI{l}) > 1, labelnamesROI{l} = ['regexp:^[' sprintf('(%s)',labelnamesROI{l}{:}) ']*$']; end
                    if orderBF > 1 && isempty(regexp(labelnamesROI{l},'(?<=bin)[ \[\(]{1,3}[1-9]','once'))
                        if ~isempty(strfind(labelnamesROI{l},'regexp')), aas_log(aap,true,'Regular expression and multiple event names do not support automatic detection of multi-order BF'); end
                        labelnamesROI{l} = ['regexp:^' strrep(labelnamesROI{l},'*','.*') 'bin [' sprintf('(%d)',1:9) ']']; % add all orders as defaults
                    end
                    if iscell(labelnamesROI{l}), labelnamesROI{l} = labelnamesROI{l}{1}; end
                end
                
                % create samples
                load(aas_getfiles_bystream(aap,'subject',subj,'roidata_firstlevel_betas'),'ROI');
                load(aas_getfiles_bystream(aap,'study',[],'valid_roi_firstlevel_betas'),'ValidROI');
                [~,roiInd] = intersect([ROI.ROIval],[ValidROI.ROIval]);
                
                betas = [ROI(roiInd).mean];
                for b = 1:size(betas,1)
                    beta = betas(b,:);
                    save(spm_file(fnSPM,'filename',sprintf('beta_%04d.mat',b)),'beta');
                end
                
                % create rois
                rois = zeros(1,size(betas,2));
                networkSpec = aas_getsetting(aap,'network.networkrois');
                for n = 1:numel(networkSpec)
                    [~,indRois] = intersect(ValidROI.ROIval,networkSpec{n});
                    rois(1,indRois) = n;
                end
                fnMaskROI = fullfile(aas_getsubjpath(aap,subj),'rois_roi.mat');
                save(fnMaskROI,'rois');
            end
            if aas_stream_has_contents(aap,'subject',subj,'connectivity') % connectivity
                labelnamesConnectivity = labelnames;
                if isstruct(labelnamesConnectivity), labelnamesConnectivity = labelnamesConnectivity.connectivity; end

                % create rois
                validROIs = readtable(aas_getfiles_bystream(aap,'study',[],'ROInames'),'ReadVariableNames',false);
                validROIs = cellfun(@(r) sscanf(r,'Atlas.cluster%d'), validROIs.Var1);
                fnCM = cellstr(aas_getfiles_bystream(aap,'subject',subj,'connectivity'));
                hdr = read_header_matlab(fnCM{1});
                dat = read_image_matlab(hdr);
                switch aas_getsetting(aap,'method')
                    case 'roi' % roi-wise each column decoded separately
                        rois = repmat([1:hdr.dim(2)],hdr.dim(1),1);
                    case 'network' % within sub-network, each sub-network decoded as whole
                        rois = zeros(hdr.dim);
                        networkSpec = aas_getsetting(aap,'network.networkrois');
                        for n = 1:numel(networkSpec)
                            [~,indRois] = intersect(validROIs,networkSpec{n});
                            rois(indRois,indRois) = n;
                        end
                end
                rois(isnan(dat)) = 0;
                fnMaskConnectivity = fullfile(aas_getsubjpath(aap,subj),'rois_connectivity.mat');
                save(fnMaskConnectivity,'rois');
            end
        end

        cfg.results.dir = fullfile(aas_getsubjpath(aap,subj), 'decoding');
        cfg.results.output = strsplit(aas_getsetting(aap,'measure'),':');
        cfg.results.overwrite = 1;
        
        %% Run
        cfg0 = cfg;
        
        % original multiclass
        if aas_stream_has_contents(aap,'subject',subj,'firstlevel_betas')
            aas_log(aap,false,'INFO: VV decoding');
            cfg = decoding_describe_data(cfg,labelnames,1:length(labelnames),regressor_names,dirSPM);
        else
            if aas_stream_has_contents(aap,'subject',subj,'roidata_firstlevel_betas') % roi
                cfgROI = decoding_describe_data_roi(cfg,labelnamesROI,1:length(labelnamesROI),regressor_names,dirSPM);
            end
            if aas_stream_has_contents(aap,'subject',subj,'connectivity') % connectivity
                cfgConnectivity = decoding_describe_data_connectivity(cfg,labelnamesConnectivity,1:length(labelnamesConnectivity), aap, subj);
            end
            if exist('cfgROI','var') && exist('cfgConnectivity','var') % connectivity + roi
                aas_log(aap,false,'INFO: joint decoding - ROI + Connectivity');
                labelnames = aas_getsetting(aap,'itemList');
                [~,indROI,indConnectivity] = intersect(ValidROI.ROIval,validROIs);
                load(fnMaskROI,'rois');
                maskROIs = rois(indROI);
                load(fnMaskConnectivity,'rois');
                maskConnectivity = rois(indConnectivity,indConnectivity);
                rois = [maskROIs; maskConnectivity];
                fnMask = fullfile(aas_getsubjpath(aap,subj),'rois.mat');
                save(fnMask,'rois');
                for f = 1:numel(cfgROI.files.name)
                    load(cfgROI.files.name{f},'beta');
                    load(cfgConnectivity.files.name{f},'sessionCM');
                    beta = [beta(indROI); sessionCM(indConnectivity,indConnectivity)];
                    save(cfgROI.files.name{f},'beta');
                end
                cfg = cfgROI;
                cfg.files.mask = fnMask;
                labelnames = labelnamesROI;
            elseif exist('cfgROI','var')
                aas_log(aap,false,'INFO: ROI decoding');
                cfg = cfgROI;
                cfg.files.mask = fnMaskROI;
                labelnames = labelnamesROI;
            elseif exist('cfgConnectivity','var')
                aas_log(aap,false,'INFO: Connectivity decoding');
                cfg = cfgConnectivity;
                cfg.files.mask = fnMaskConnectivity;
                labelnames = labelnamesConnectivity;
            end
        end
        aap = aas_desc_outputs(aap,'subject',subj,'samples',cfg.files.name);
        cfg0.files.mask = cfg.files.mask;

        cfg.design = make_design_cv(cfg);
        decoding(cfg);
        
        dopairwise = aas_getsetting(aap,'dopairwise') && numel(labelnames) > 2;
        doonevsall = aas_getsetting(aap,'doonevsall') && numel(labelnames) > 2;
        
        % pairwise
        if dopairwise
            % - create labelcombinations
            labelcmbsel = [reshape(repmat(1:numel(labelnames),numel(labelnames),1),1,[]); repmat(1:numel(labelnames),1,numel(labelnames))]';
            labelcmbsel(diff(labelcmbsel,[],2) == 0,:) = [];
            labelcmbsel = unique(cell2mat(arrayfun(@(c) sort(labelcmbsel(c,:))', 1:size(labelcmbsel,1), 'UniformOutput', false))','rows');
            
            % - run pairwise decodings
            for c = 1:size(labelcmbsel,1)
                cfg = cfg0;
                cfg.results.filestart = sprintf('%spw%02d',cfg.results.filestart,c); % assume no more than 99 pairwise comparisons, i.e. 14 labels (see also lines 255 and 263)
                if aas_stream_has_contents(aap,'subject',subj,'firstlevel_betas')
                    cfg = decoding_describe_data(cfg,labelnames(labelcmbsel(c,:)),1:2,regressor_names,dirSPM);
                elseif aas_stream_has_contents(aap,'subject',subj,'roidata_firstlevel_betas') % roi
                    cfg = decoding_describe_data_roi(cfg,labelnames(labelcmbsel(c,:)),1:2,regressor_names,dirSPM);
                elseif aas_stream_has_contents(aap,'subject',subj,'connectivity') % connectivity 
                    cfg = decoding_describe_data_connectivity(cfg,labelnames(labelcmbsel(c,:)),1:2, aap, subj);
                end
                cfg.design = make_design_cv(cfg);
                decoding(cfg);
            end
            aap = aas_desc_outputs(aap,'subject',subj,'settings_pairwise',spm_select('FPList',cfg0.results.dir,['^' cfg0.results.filestart 'pw[0-9]{2}_cfg.mat']));
        end
        
        % one-vs-all decodings
        if doonevsall
            for c = 1:numel(labelnames)
                labels = ones(1,numel(labelnames));
                labels(c) = 2;
                cfg = cfg0;
                cfg.results.filestart = sprintf('%sova%02d',cfg.results.filestart,c); % assume no more than 99 pairwise comparisons, i.e. 14 labels (see also lines 255 and 263)
                if aas_stream_has_contents(aap,'subject',subj,'firstlevel_betas')
                    cfg = decoding_describe_data(cfg,labelnames,labels,regressor_names,dirSPM);
                elseif aas_stream_has_contents(aap,'subject',subj,'roidata_firstlevel_betas') % roi
                    cfg = decoding_describe_data_roi(cfg,labelnames,labels,regressor_names,dirSPM);
                elseif aas_stream_has_contents(aap,'subject',subj,'connectivity') % connectivity 
                    cfg = decoding_describe_data_connectivity(cfg,labelnames,labels, aap, subj);
                end
                cfg.design = make_design_cv(cfg);
                cfg.design.unbalanced_data = 'ok';
                decoding(cfg);
            end
            aap = aas_desc_outputs(aap,'subject',subj,'settings_onevsall',spm_select('FPList',cfg0.results.dir,['^' cfg0.results.filestart 'ova[0-9]{2}_cfg.mat']));
        end
        
        aap = aas_desc_outputs(aap,'subject',subj,'settings',fullfile(cfg0.results.dir,[cfg0.results.filestart '_cfg.mat']));
        aap = aas_desc_outputs(aap,'subject',subj,'mask',cfg0.files.mask);
                
        for o = 1:numel(cfg0.results.output)
            aap = aas_desc_outputs(aap,'subject',subj,cfg0.results.output{o},spm_select('FPList',cfg0.results.dir,['^' cfg0.results.filestart '_' cfg0.results.output{o} '[_a-z]*\.[(mat)(nii)]{1}']));
            if dopairwise
                aap = aas_desc_outputs(aap,'subject',subj,[cfg0.results.output{o} '_pairwise'],spm_select('FPList',cfg0.results.dir,['^' cfg0.results.filestart 'pw[0-9]{2}_' cfg0.results.output{o} '[_a-z]*\.[(mat)(nii)]{1}']));
            end
            if doonevsall
                aap = aas_desc_outputs(aap,'subject',subj,[cfg0.results.output{o} '_onevsall'],spm_select('FPList',cfg0.results.dir,['^' cfg0.results.filestart 'ova[0-9]{2}_' cfg0.results.output{o} '[_a-z]*\.[(mat)(nii)]{1}']));
            end
        end
        
        %% Cleanup
        TDT.unload;
        
    case 'checkrequirements'
        if subj == 1 % only once
            if ~aas_cache_get(aap,'tdt'), aas_log(aap,true,'TDT is not found'); end
            
            dopairwise = aas_getsetting(aap,'dopairwise') && numel(aas_getsetting(aap,'itemList')) > 2;
            doonevsall = aas_getsetting(aap,'doonevsall') && numel(aas_getsetting(aap,'itemList')) > 2;
            if dopairwise, aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'append','settings_pairwise','output'); end
            if doonevsall, aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'append','settings_onevsall','output'); end
            
            out = aas_getstreams(aap,'output'); out = setdiff(out,{'settings' 'settings_pairwise' 'samples' 'mask'});
            meas = aas_getsetting(aap,'measure'); meas = strsplit(meas,':');
            if isempty(meas{1}), aas_log(aap,true,'no measure specified'); end
            for s = 1:numel(meas)
                if strcmp(out{s},meas{s}), continue; end
                if s == 1
                    aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'result',meas{s},'output');
                else
                    aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'append',meas{s},'output');
                end
                if dopairwise, aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'append',[meas{s} '_pairwise'],'output'); end
                if doonevsall, aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'append',[meas{s} '_onevsall'],'output'); end
                aas_log(aap,false,['INFO: ' aap.tasklist.currenttask.name ' output stream: ''' meas{s} '''']);
            end
        end
end
end

function cfg = decoding_describe_data_roi(cfg,labelnames,labels, regressor_names, dirSPM)
cfg.files.name = {};
cfg.files.chunk = [];
cfg.files.label = [];
cfg.files.labelname = {};
cfg.files.set = [];
cfg.files.xclass = [];
cfg.files.descr = {};
fns = cellstr(spm_select('FPList',dirSPM,'^beta_[0-9]{4}.mat'));
for l = 1:numel(labelnames)
    fnInd = cellfun(@(r) ~isempty(regexp(r,wildcard2regexp(labelnames{l}), 'once')), regressor_names(1,:));
    fnsLabel = fns(fnInd);
    
    cfg.files.name = vertcat(cfg.files.name,fnsLabel);
    cfg.files.chunk = vertcat(cfg.files.chunk,...
        [regressor_names{2,fnInd}]');
    cfg.files.label = vertcat(cfg.files.label,repmat(labels(l),numel(fnsLabel),1));
    cfg.files.labelname = vertcat(cfg.files.labelname,...
        repmat(labelnames(l),numel(fnsLabel),1));
    cfg.files.descr = vertcat(cfg.files.descr,regressor_names(3,fnInd)');
end
end

function cfg = decoding_describe_data_connectivity(cfg,labelnames,labels, aap, subj)

cfg.files.name = {};
cfg.files.chunk = [];
cfg.files.label = [];
cfg.files.labelname = {};
cfg.files.set = [];
cfg.files.xclass = [];
cfg.files.descr = {};
for sess = aap.acq_details.selected_sessions
    fns = cellstr(aas_getfiles_bystream(aap,'session',[subj sess],'connectivity'));    
    for l = 1:numel(labelnames)
        fnptrnLabel = spm_file(labelnames{l}{1},'ext','mat');
        fnsLabel = cellstr(get_filenames(cfg.software,spm_file(fns{1},'path'),fnptrnLabel));
        if isempty(fnsLabel{1}), fnsLabel = cellstr(get_filenames(cfg.software,spm_file(fns{1},'path'),['*' fnptrnLabel])); end
        
        cfg.files.name = vertcat(cfg.files.name,fnsLabel);
        cfg.files.chunk = vertcat(cfg.files.chunk,...
            sess*ones(numel(fnsLabel),1));
        cfg.files.label = vertcat(cfg.files.label,labels(l)*ones(numel(fnsLabel),1));
        cfg.files.labelname = vertcat(cfg.files.labelname,...
            repmat(cellstr(strjoin(labelnames{l},'+')),numel(fnsLabel),1));
        cfg.files.descr = vertcat(cfg.files.descr,spm_file(fns(l),'basename'));
    end
end

% reorder to match with the others (labels -> chunks)
[~,ind] = sortrows([cfg.files.label cfg.files.chunk]);
cfg.files.name = cfg.files.name(ind);
cfg.files.chunk = cfg.files.chunk(ind);
cfg.files.label = cfg.files.label(ind);
cfg.files.labelname = cfg.files.labelname(ind);
cfg.files.descr = cfg.files.descr(ind);
end