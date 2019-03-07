function [aap,resp]=aamod_CoSMoMVPA(aap,task,subj)
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
        %% Prepare data
        TASKS = textscan(aas_getsetting(aap,'tasks'),'%s','Delimiter',':'); TASKS = TASKS{1}';
        
        RSAROOT = fullfile(aas_getsubjpath(aap,subj),'RSA');
        aas_makedir(aap,RSAROOT);
        
        inps = aas_getstreams(aap,'input');
        inps = inps(logical(cellfun(@(x) exist(aas_getinputstreamfilename(aap,'subject',subj,x),'file'), inps)));
        struct_fn = aas_getfiles_bystream(aap,'subject',subj,inps{1});
        
        fnMask = cellfun(@(x) aas_getfiles_bystream(aap,'subject',subj,x), inps(cell_index(inps,'mask')),'UniformOutput',false);
        fnSPM = cellfun(@(x) aas_getfiles_bystream(aap,'subject',subj,x), inps(cellfun(@(x) ~isempty(regexp(x,'firstlevel_spm$', 'once')), inps)),'UniformOutput',false);
        fnTmaps = cellfun(@(x) aas_getfiles_bystream(aap,'subject',subj,x), inps(cell_index(inps,'firstlevel_spmts')),'UniformOutput',false);
        
        if numel(fnMask) > 1
            brain_mask = spm_imcalc(spm_vol(char(fnMask)),fullfile(RSAROOT,'brain_mask.nii'),'min(X)',{1});
        else
            brain_mask.fname = char(fnMask);
        end
        
        ITEMS = aas_getsetting(aap,'itemList');
        fnT = {};
        for run = 1:numel(fnSPM)
            load(fnSPM{run});
            conNames = {SPM.xCon.name}';
            for it = 1:numel(ITEMS)
                ITEMS{it} = cellstr(ITEMS{it});
                fnT{end+1} = fnTmaps{run}(sum(cellfun(@(x) cell_index(conNames,x), ITEMS{it})),:);
            end
        end
        spm_file_merge(fnT,fullfile(RSAROOT,'glm_T_stats_perrun.nii'));
        
        %% Initialise Cosmo
        oldPath = path;
        cosmo_set_path
        cosmo_check_external('-tic');
        
        if cell_index(TASKS, 'RSA') && isempty(aas_getsetting(aap,'bsMatrix')), TASKS(cell_index(TASKS, 'RSA')) = []; end
        
        for t = TASKS
            % Data
            switch t{1}
                case 'RSA'
                    ds=cosmo_fmri_dataset(fullfile(RSAROOT,'glm_T_stats_perrun.nii'),'mask',brain_mask.fname,...
                        'targets',repmat(1:numel(ITEMS),1,numel(fnSPM))');
                    ds=cosmo_fx(ds, @(x)mean(x,1), 'targets', 1);
                    ds.sa.labels=cellfun(@(x) x{1}, ITEMS, 'UniformOutput', false)';
                    ds.sa.set=(1:numel(ITEMS))';
                case 'C'
                    ds=cosmo_fmri_dataset(fullfile(RSAROOT,'glm_T_stats_perrun.nii'),'mask',brain_mask.fname,...
                        'targets',repmat(1:numel(ITEMS),1,numel(fnSPM))','chunks',floor(((1:numel(ITEMS)*numel(fnSPM))-1)/numel(ITEMS))+1);
                    ds.sa.labels=cellfun(@(x) x{1}, repmat(ITEMS,1,numel(fnSPM)), 'UniformOutput', false)';
                    ds.sa.set=repmat((1:numel(ITEMS))',numel(fnSPM),1);
            end            
            
            % Data            
            cosmo_check_dataset(ds);
            
            % Searchlight
            nbrhood=cosmo_spherical_neighborhood(ds,'count',aas_getsetting(aap,'searchlightVox'));
            
            % Model
            target_dsm = 'not specified';
            switch t{1}
                case 'RSA'
                    target_dsm=importdata(aas_getsetting(aap,'bsMatrix'));
                    measure=@cosmo_target_dsm_corr_measure;
                    measure_args = aas_getsetting(aap,'RSAsettings');
                    measure_args.target_dsm=target_dsm;
                case 'C'
                    measure_args = aas_getsetting(aap,'Csettings');
                    measure_args.classifier = str2func(['cosmo_classify_' lower(measure_args.classifier)]);
                    measure_args.partitions = cosmo_balance_partitions(cosmo_nfold_partitioner(ds), ds);
                    measure=@cosmo_crossvalidation_measure;
            end  
            
            %% Info
            aas_log(aap,false,'INFO:Dataset input:'); cosmo_disp(ds);
            aas_log(aap,false,'INFO:Searchlight neighborhood definition:'); cosmo_disp(nbrhood);
            aas_log(aap,false,'INFO:Target DSM:'); disp(target_dsm);
            
            %         imagesc(target_dsm)
            %         set(gca,'XTick',1:size(ds.samples,1),'XTickLabel',ds.sa.labels,...
            %             'YTick',1:size(ds.samples,1),'YTickLabel',ds.sa.labels)
            
            %% Run
            ds_rsm_behav=cosmo_searchlight(ds,nbrhood,measure,measure_args);
            
            %         cosmo_plot_slices(ds_rsm_behav);
            
            % store results
            rsa_fn=fullfile(RSAROOT,[t{1} 'map.nii']);
            cosmo_map2fmri(ds_rsm_behav,rsa_fn);
            
            %% Cleanup
            path(oldPath);
            
            aap=aas_desc_outputs(aap,'subject',subj,[t{1} 'map'],rsa_fn);
        end
    case 'checkrequirements'
        
    otherwise
        if isempty(which('cosmo_set_path')), aas_log(aap,true,sprintf('CoSMoMVPA cannot be found!\n Make sure you add <CoSMoMVPA directory>/mvpa to aap.directory_conventions.spmtoolsdir.')); end
        reqInps = {'firstlevel_brainmask' 'firstlevel_spm' 'firstlevel_spmts'};
        inps = aas_getstreams(aap,'input');
        missingInps = reqInps(cellfun(@(x) ~any(cell_index(inps,x)), reqInps));
        if ~isempty(missingInps), aas_log(aap,true,['Inputs not specified:' sprintf(' %s',missingInps{:})]); end
end