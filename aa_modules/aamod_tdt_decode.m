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
        cfg = TDT.defaults;
        cfg.analysis = aas_getsetting(aap,'method');
        
        % Disable plots
        cfg.plot_selected_voxels = 0;
        cfg.plot_design = 0;
        
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
        
        fnMask = cellfun(@(x) aas_getfiles_bystream(aap,'subject',subj,x), inps(contains(inps,{'mask','rois'})),'UniformOutput',false);
        fnSPM = aas_getfiles_bystream(aap,'subject',subj,'firstlevel_spm');
        dat = load(fnSPM);
        dirSPM = fullfile(aas_getsubjpath(aap,subj),aap.directory_conventions.stats_singlesubj);
        orderBF = dat.SPM.xBF.order;
        
        % check mask
        % ensure overlap with data
        fnBeta = cellstr(aas_getfiles_bystream(aap,'subject',subj,'firstlevel_betas'));
        Ydata = spm_read_vols(spm_vol(fnBeta{1}));
        if (numel(fnMask) > 1) && ~strcmp(cfg.analysis,'roi')
            brain_mask = spm_imcalc(spm_vol(char(fnMask)),fullfile(TASKROOT,'brain_mask.nii'),'min(X)',{1});
            fnMask = {brain_mask.fname};
        end
        for f = 1:numel(fnMask)
            V = spm_vol(fnMask{f});
            Y = spm_read_vols(V);
            mask_num = unique(Y(:));
            mask_num(isnan(mask_num)|mask_num==0) = [];
            if numel(mask_num) > 1 && ~strcmp(cfg.analysis,'roi')
                fnMask{f} = spm_file(V.fname,'prefix','b_');
                Y = Y > 0.5;
            end
            fnMask{f} = spm_file(fnMask{f},'prefix','ok_');
            Y(isnan(Ydata)) = 0;
            V.fname = fnMask{f};
            V.pinfo = [1 0 0]';
            spm_write_vol(V,Y);
        end
        
        cfg.files.mask = fnMask;
        
        switch cfg.analysis
            case 'searchlight'
                cfg.searchlight = aas_getsetting(aap,'searchlight');
        end
        
        labelnames = aas_getsetting(aap,'itemList');
        for l = 1:numel(labelnames)
            if numel(labelnames{l}) > 1, labelnames{l} = ['regexp:^[' sprintf('(%s)',labelnames{l}{:}) ']*$']; end
            if orderBF > 1 && isempty(regexp(labelnames{l},'(?<=bin)[ \[\(]{1,3}[1-9]','once'))
                if ~isempty(strfind(labelnames{l},'regexp')), aas_log(aap,true,'Regular expression and multiple event names do not support automatic detection of multi-order BF'); end
                labelnames{l} = ['regexp:^' strrep(labelnames{l},'*','.*') 'bin [' sprintf('(%d)',1:9) ']']; % add all orders as defaults
            end
            if iscell(labelnames{l}), labelnames{l} = labelnames{l}{1}; end
        end
 
        cfg.results.dir = fullfile(aas_getsubjpath(aap,subj), 'decoding');
        cfg.results.output = strsplit(aas_getsetting(aap,'measure'),':');
        
        %% Run
        regressor_names = design_from_spm(dirSPM);
        
        cfg0 = cfg;
        if aas_getsetting(aap,'pairwise')
            % create labelcombinations
            labelcmbsel = [reshape(repmat(1:numel(labelnames),numel(labelnames),1),1,[]); repmat(1:numel(labelnames),1,numel(labelnames))]';
            labelcmbsel(diff(labelcmbsel,[],2) == 0,:) = []; 
            labelcmbsel = unique(cell2mat(arrayfun(@(c) sort(labelcmbsel(c,:))', 1:size(labelcmbsel,1), 'UniformOutput', false))','rows');
            
            % run pairwise decodings
            for c = 1:size(labelcmbsel,1)
                cfg = cfg0;
                cfg.results.filestart = sprintf('%s%02d',cfg.results.filestart,c);
                cfg = decoding_describe_data(cfg,labelnames(labelcmbsel(c,:)),1:2,regressor_names,dirSPM);
                cfg.design = make_design_cv(cfg);
                decoding(cfg);
            end
        else
            cfg = decoding_describe_data(cfg,labelnames,1:length(labelnames),regressor_names,dirSPM);
            cfg.design = make_design_cv(cfg);
            decoding(cfg);
        end
        
        aap = aas_desc_outputs(aap,'subject',subj,'settings',spm_select('FPList',cfg0.results.dir,[cfg0.results.filestart '.*_cfg.mat']));
        aap = aas_desc_outputs(aap,'subject',subj,'mask',cfg0.files.mask);
                
        for o = 1:numel(cfg0.results.output)
            aap = aas_desc_outputs(aap,'subject',subj,cfg0.results.output{o},spm_select('FPList',cfg0.results.dir,['^' cfg0.results.filestart '[0-9]*_' cfg0.results.output{o} '[_a-z]*\.[(mat)(nii)]{1}']));
        end
        
        %% Cleanup
        TDT.unload;
        
    case 'checkrequirements'
        if ~aas_cache_get(aap,'tdt'), aas_log(aap,true,'TDT is not found'); end
                
        out = aas_getstreams(aap,'output'); out(1:2) = []; % settings and mask
        meas = aas_getsetting(aap,'measure'); meas = strsplit(meas,':');
        if isempty(meas{1}), aas_log(aap,true,'no measure specified'); end
        for s = 1:numel(meas)
            if strcmp(out{s},meas{s}), continue; end
            if s == 1
                aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'result',meas{s},'output');
            else
                aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'append',meas{s},'output');
            end
            aas_log(aap,false,['INFO: ' aap.tasklist.currenttask.name ' output stream: ''' meas{s} '''']);
        end
end
end