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
        cfg.analysis = 'searchlight';
        
        % Disable plots
        cfg.plot_selected_voxels = 0;
        cfg.plot_design = 0;
        
        %% Prepare
        inps = aas_getstreams(aap,'input');
        inps = inps(logical(cellfun(@(x) exist(aas_getinputstreamfilename(aap,'subject',subj,x),'file'), inps)));
        
        fnMask = cellfun(@(x) aas_getfiles_bystream(aap,'subject',subj,x), inps(cell_index(inps,'mask')),'UniformOutput',false);
        fnSPM = aas_getfiles_bystream(aap,'subject',subj,'firstlevel_spm');
        dat = load(fnSPM);
        dirSPM = fullfile(aas_getsubjpath(aap,subj),aap.directory_conventions.stats_singlesubj);
        orderBF = dat.SPM.xBF.order;
        
        % check mask
        if numel(fnMask) > 1
            brain_mask = spm_imcalc(spm_vol(char(fnMask)),fullfile(TASKROOT,'brain_mask.nii'),'min(X)',{1});
            fnMask = brain_mask.fname;
        else
            fnMask = char(fnMask);
        end
        V = spm_vol(fnMask);
        Y = spm_read_vols(V); 
        mask_num = unique(Y(:));
        mask_num(isnan(mask_num)|mask_num==0) = []; 
        if any(mask_num~=1)
            fnMask = spm_file(V.fname,'prefix','b');
            V.fname = fnMask;
            V.pinfo = [1 0 0]';
            Y = Y > 0.5;
            spm_write_vol(V,Y);
        end
        
        cfg.files.mask = fnMask;
        
        cfg.searchlight = aas_getsetting(aap,'searchlight');
        
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
        
        cfg = decoding_describe_data(cfg,labelnames,1:length(labelnames),regressor_names,dirSPM);
        
        cfg.design = make_design_cv(cfg);

        [~, cfg] = decoding(cfg);
        
        aap = aas_desc_outputs(aap,'subject',subj,'settings',fullfile(cfg.results.dir,[cfg.results.filestart '_cfg.mat']));
        aap = aas_desc_outputs(aap,'subject',subj,'mask',cfg.files.mask);
                
        for o = 1:numel(cfg.results.output)
            aap = aas_desc_outputs(aap,'subject',subj,cfg.results.output{o},spm_select('FPList',cfg.results.dir,['^' cfg.results.resultsname{o} '\.[(mat)(nii)]{1}']));
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