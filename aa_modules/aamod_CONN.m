function [aap,resp]=aamod_CONN_validate(aap,task,subj)
resp='';

switch task
%     case 'report'
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
        [~, CONN] = aas_cache_get(aap, 'conn');
        CONN.load;
        global CONN_x
        
        roiNames = {'GM' 'WM' 'CSF' 'Atlas'};
        
        CONNROOT = fullfile(aas_getsubjpath(aap,subj), 'conn');
        aas_makedir(aap,CONNROOT);
        
        % Setup
        bcpvar = aap.options.verbose;
        aap.options.verbose = -1;
        fnAnat = aas_getfiles_bystream(aap,'subject',subj,'structural');
        for sess = 1:aas_getN_bydomain(aap,'session',subj)
            fnFunc{sess} = aas_getfiles_bystream(aap,'session',[subj sess],'epi');
        end        
        fnROIs.(roiNames{1}) = aas_getfiles_bystream(aap,'subject',subj,'native_grey_mask');
        fnROIs.(roiNames{2}) = aas_getfiles_bystream(aap,'subject',subj,'native_white_mask');
        fnROIs.(roiNames{3}) = aas_getfiles_bystream(aap,'subject',subj,'native_csf_mask');
        fnROIs.(roiNames{4}) = aas_getfiles_bystream(aap,'subject',subj,'rois');
        fnSPM = aas_getfiles_bystream(aap,'subject',subj,'firstlevel_spm');
        aap.options.verbose = bcpvar;
        analyses = aas_getsetting(aap,'analysis');
        nSess = aas_getN_bydomain(aap,'session',[]);
        
        % Heuristics
        if ~isempty(fnROIs.Atlas)
            Y = spm_read_vols(spm_vol(fnROIs.Atlas));
            multiROI = sum(unique(Y)~=0)>1;
            ROIonly = all(cellfun(@isempty, {analyses.roival}));
            doPPI = ~cellfun(@isempty, {analyses.condition});
        else
            multiROI = false;
            ROIonly = false;
            doPPI = false;
        end
        if ROIonly, processNameSfx = '_roi'; end
        
        % Initialize
        CONN_x.Setup.analysisunits = 2; % 1 = PSC, 2 = raw
        CONN_x.Setup.outputfiles = [0 0 0 0 0 0];
        
        CONN_x.filename = fullfile(CONNROOT);
        CONN_x.folders.data = fullfile(CONN_x.filename,'data');
        CONN_x.folders.preprocessing = fullfile(CONN_x.filename,'preprocessing');
        CONN_x.folders.qa = fullfile(CONN_x.filename,'preprocessing','qa');
        CONN_x.folders.bookmarks = fullfile(CONN_x.filename,'results','bookmarks');
        CONN_x.folders.firstlevel = fullfile(CONN_x.filename,'results','firstlevel');
        CONN_x.folders.firstlevel_vv = fullfile(CONN_x.filename,'results','firstlevel_vv');
        CONN_x.folders.firstlevel_dyn = fullfile(CONN_x.filename,'results','firstlevel_dyn');
        CONN_x.folders.secondlevel = fullfile(CONN_x.filename,'results','secondlevel');
        aas_makedir(aap,CONN_x.folders.data);
        aas_makedir(aap,CONN_x.folders.preprocessing);
        
        % Data
        if ~isempty(fnAnat)
            aas_log(aap,false,'INFO: structural found');
            CONN_x.Setup.structural{1} = repmat({conn_file(fnAnat)},1,nSess);
            CONN_x.Setup.structural_sessionspecific = 1;
        end
        CONN_x.Setup.functional = {cellfun(@conn_file, fnFunc, 'UniformOutput', false)};
        CONN_x.Setup.rois.files{1} = {};
        for r = roiNames
            if ~isempty(fnROIs.(r{1}))
                aas_log(aap,false,['INFO: ' r{1} ' mask found']);
                CONN_x.Setup.rois.names(end+1) = r; 
                CONN_x.Setup.rois.files{1}{end+1} = repmat({conn_file(fnROIs.(r{1}))},1,nSess);                
            end
        end
        CONN_x.Setup.rois.names(end+1) = {''};
        nROIs = numel(CONN_x.Setup.rois.names)-1;
        CONN_x.Setup.rois.dimensions = repmat({1},1,nROIs);
        [~, indTissue] = intersect(CONN_x.Setup.rois.names,{'WM' 'CSF'});
        CONN_x.Setup.rois.dimensions(indTissue) = {16};
        CONN_x.Setup.rois.mask = false(1,nROIs);
        CONN_x.Setup.rois.subjectspecific = true(1,nROIs);
        CONN_x.Setup.rois.sessionspecific = true(1,nROIs);
        CONN_x.Setup.rois.multiplelabels = false(1,nROIs);
        CONN_x.Setup.rois.multiplelabels(strcmp(CONN_x.Setup.rois.names,'Atlas')) = multiROI;
        CONN_x.Setup.rois.regresscovariates = false(1,nROIs);
        CONN_x.Setup.rois.regresscovariates(indTissue) = true;
        CONN_x.Setup.rois.unsmoothedvolumes = true(1,nROIs);
        CONN_x.Setup.rois.weighted = false(1,nROIs);        
        if ~isempty(fnSPM)
            aas_log(aap,false,'INFO: SPM found');
            conn_importspm(fnSPM,...
                'addfunctional',false,...
                'addconditions',true,...
                'breakconditionsbysession',true,...
                'addrestcondition',false,...
                'keeppreviousconditions',false,...
                'addcovariates',true,...
                'addrealignment',false,...
                'addartfiles',false);
        end        
        conn_process(['setup' processNameSfx]);
        
        % Denoising
        CONN_x.filename = fullfile(CONNROOT);
        [~,indConf] = intersect(CONN_x.Preproc.variables.names, aas_getsetting(aap,'denoising.confounds')); indConf = sort(indConf);
        CONN_x.Preproc.confounds.names = CONN_x.Preproc.variables.names(indConf);
        CONN_x.Preproc.confounds.filter = repmat({0},1,numel(CONN_x.Preproc.confounds.names));
        CONN_x.Preproc.confounds.types = {};
        CONN_x.Preproc.confounds.power = {};
        CONN_x.Preproc.confounds.deriv = {};
        CONN_x.Preproc.confounds.dimensions = {};
        
        CONN_x.Preproc.filter = aas_getsetting(aap,'denoising.filter');
        if isnan(CONN_x.Preproc.filter(1))
            CONN_x.Preproc.filter(1) = 1/aas_getsetting(aas_setcurrenttask(aap,aas_getsourcestage(aap,'aamod_firstlevel_model')),'highpassfilter');
        end
        CONN_x.Preproc.regbp = 2; % 1 = filter, then regress, 2 = Simultaneous regression and filtering
        CONN_x.Preproc.detrending = 1;        
        conn_process(['denoising' processNameSfx]);        
        
        % Setup analyses
        for a = 1:numel(analyses)
            CONN_x.Analyses(a).name = analyses(a).name;
            CONN_x.Analyses(a).conditions = intersect(CONN_x.Setup.conditions.names,analyses(a).condition);
            if ~isempty(analyses(a).roival)
                CONN_x.Analyses(a).regressors.names = intersect(CONN_x.Analysis_variables.names,arrayfun(@(r) sprintf('%s.cluster%03d','Atlas',r), analyses(a).roival, 'UniformOutput', false));
            else
                CONN_x.Analyses(a).regressors.names = CONN_x.Analysis_variables.names(contains(CONN_x.Analysis_variables.names,'Atlas'));
            end
            CONN_x.Analyses(a).regressors.dimensions = {};
            CONN_x.Analyses(a).regressors.deriv = {};
            CONN_x.Analyses(a).regressors.fbands = {};
        end
        conn_process('analyses_seedsetup');
        
        % Run analyses
        CONN_x.gui = [];
        for a = 1:numel(analyses)
            if numel(CONN_x.Analyses(a).regressors.names) == 1 % 'Seed-to-Voxel'
                CONN_x.Analyses(a).type = 2;
            else % 'ROI-to-ROI'
                CONN_x.Analyses(a).type = 1;
            end
            CONN_x.Analyses(a).measure = find(strcmp(strsplit(aap.schema.tasksettings.aamod_CONN(1).analysis.measure.ATTRIBUTE.options,'|'),analyses(a).measure));
            CONN_x.Analyses(a).modulation = doPPI(a);
            CONN_x.Analyses(a).weight = find(strcmp(strsplit(aap.schema.tasksettings.aamod_CONN(1).analysis.weight.ATTRIBUTE.options,'|'),analyses(a).weight));
            [~,indCond] = intersect(CONN_x.Setup.conditions.allnames,CONN_x.Analyses(a).conditions);
            CONN_x.gui.conditions = indCond';
            conn_process('analyses_gui_seedandroi',a)
        end
        
        save(CONNROOT,'CONN_x');
        aap = aas_desc_outputs(aap,'subject',subj,'settings',[CONNROOT '.mat']);
        
        % Results/summary
        fnOutSess = repmat({{}},1,nSess); fnOutSubj = {};
        for a = 1:numel(CONN_x.Analyses)
            CONN_x.Results(1).name = CONN_x.Analyses(a).name;
            CONN_x.Results.foldername = fullfile(CONN_x.Results.name);
            CONN_x.Results.display = 0;
            CONN_x.Results.xX.nsubjecteffects = 1;
            CONN_x.Results.xX.csubjecteffects = 1;
            CONN_x.Results.xX.nsubjecteffectsbyname = {'AllSubjects'};
            [~, indCond] = intersect(CONN_x.Setup.conditions.allnames,CONN_x.Analyses(a).conditions);
            CONN_x.Results.xX.nconditions = indCond;
            CONN_x.Results.xX.cconditions = ones(1,numel(CONN_x.Analyses(a).conditions));
            CONN_x.Results.xX.nconditionsbyname = CONN_x.Analyses(a).conditions;
%             CONN_x.Results.xX.modeltype = 2; % 1 = RFX, 2 = FFX - not working!
            CONN_x.Analysis = a;
            conn_process(['results' processNameSfx]);
            
            % Session
            for c = 1:numel(CONN_x.Results.xX.nconditions)
                sessionCM = zeros(numel(CONN_x.Analyses(a).sources),numel(CONN_x.Analyses(a).sources),1);
                fn = fullfile(CONN_x.folders.firstlevel,CONN_x.Analyses(a).name,sprintf('resultsROI_Condition%03d.mat',CONN_x.Results.xX.nconditions(c)));
                dat = load(fn);
                sessionCM(:,:,1) = dat.Z;
                sessN = str2double(regexp(CONN_x.Results.xX.nconditionsbyname{c},'(?<=Session)[0-9]','match'));
                fnOutSess{sessN}(end+1) = cellstr(spm_file(fullfile(aas_getsesspath(aap,subj,sessN),CONN_x.Analyses(a).name),'ext','mat'));
                save(fnOutSess{sessN}{end},'sessionCM');
            end
            % Subject
            dat = load(fullfile(CONN_x.folders.secondlevel,CONN_x.Analyses(a).name,'ROI.mat'));
            subjectCM = vertcat(dat.ROI.h);
            fnOutSubj(end+1) = cellstr(spm_file(fullfile(aas_getsubjpath(aap,subj),CONN_x.Analyses(a).name),'ext','mat'));
            save(fnOutSubj{end},'subjectCM');            
        end
        for sess = aap.acq_details.selected_sessions
            aap = aas_desc_outputs(aap,'session',[subj sess],'connectivity',fnOutSess{sess});
        end
        aap = aas_desc_outputs(aap,'subject',subj,'connectivity',fnOutSubj);   
        
        CONN.unload;
        
    case 'checkrequirements'
        if ~aas_cache_get(aap,'conn'), aas_log(aap,true,'CONN is not found'); end
end
