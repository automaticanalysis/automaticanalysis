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
        
        TASKROOT = fullfile(aas_getsubjpath(aap,subj),'CoSMo');
        aas_makedir(aap,TASKROOT);
        
        inps = aas_getstreams(aap,'input');
        inps = inps(logical(cellfun(@(x) exist(aas_getinputstreamfilename(aap,'subject',subj,x),'file'), inps)));
        %struct_fn = aas_getfiles_bystream(aap,'subject',subj,inps{1});
        
        fnMask = cellfun(@(x) aas_getfiles_bystream(aap,'subject',subj,x), inps(cell_index(inps,'mask')),'UniformOutput',false);
        fnSPM = aas_getfiles_bystream(aap,'subject',subj,'firstlevel_spm');
        fnTmaps = cellstr(aas_getfiles_bystream(aap,'subject',subj,'firstlevel_spmts'));
        
        if numel(fnMask) > 1
            brain_mask = spm_imcalc(spm_vol(char(fnMask)),fullfile(TASKROOT,'brain_mask.nii'),'min(X)',{1});
            fnMask = brain_mask.fname;
        else
            fnMask = char(fnMask);
        end
        dat = load(fnSPM);
        firstlevelSPM = dat.SPM;
        
        ITEMS = cellfun(@cellstr, aas_getsetting(aap,'itemList'),'UniformOutput', false);
        indTmaps = [];
        inpdep = aap.internal.inputstreamsources{aap.tasklist.currenttask.modulenumber}.stream;
        firstlevelaap = aas_setcurrenttask(aap, inpdep(strcmp({inpdep.name},'firstlevel_spm')).sourcenumber);
        consubj = aas_getsetting(firstlevelaap,'contrasts',find(arrayfun(@(x) strcmp(x.subject,aas_getsubjname(aap,subj)),aas_getsetting(firstlevelaap,'contrasts'))));
        consubj = consubj.con;
        for sess = aap.acq_details.selected_sessions
            consess = consubj(arrayfun(@(x) strcmp(x.session.names{1},aap.acq_details.sessions(sess).name),consubj));
            consessnames = cellfun(@(x) strrep(x,['_' aap.acq_details.sessions(sess).name],''),{consess.name},'UniformOutput', false);
            cons = consess(cellfun(@(i) cell2mat(cellfun(@(ii) find(strcmp(consessnames,ii)), i,'UniformOutput', false)), ITEMS));
            indTmaps = [indTmaps cellfun(@(x) find(strcmp({firstlevelSPM.xCon.name},x)), {cons.name})];
        end
        spm_file_merge(fnTmaps(indTmaps),fullfile(TASKROOT,'glm_T_stats_perrun.nii'));
        fnTstats = fullfile(TASKROOT,'glm_T_stats_perrun.nii');
        nSess = numel(aap.acq_details.selected_sessions);
        
        %% Initialise Cosmo
        [junk, MVPA] = aas_cache_get(aap,'cosmomvpa');
        MVPA.load;
        
        if cell_index(TASKS, 'RSA') && isempty(aas_getsetting(aap,'bsMatrix')), TASKS(cell_index(TASKS, 'RSA')) = []; end
        
        for t = TASKS
            % Data
            switch t{1}
                case 'RSA'
                    try
                        ds=cosmo_fmri_dataset(fnTstats,'mask',fnMask,...
                            'targets',repmat(1:numel(ITEMS),1,nSess)');
                    catch E
                        if ~isempty(strfind(E.message,'sform and qform in the NIfTI header differ'))
                            ds=cosmo_fmri_dataset(fix_qform(fnTstats),'mask',fix_qform(fnMask),...
                                'targets',repmat(1:numel(ITEMS),1,nSess)');
                        else
                            aas_log(aap,true,E.message)
                        end
                    end
                    ds=cosmo_fx(ds, @(x)mean(x,1), 'targets', 1);
                    ds.sa.labels=cellfun(@(x) x{1}, ITEMS, 'UniformOutput', false)';
                    ds.sa.set=(1:numel(ITEMS))';
                case 'C'
                    try
                        ds=cosmo_fmri_dataset(fnTstats,'mask',fnMask,...
                            'targets',repmat(1:numel(ITEMS),1,nSess)','chunks',floor(((1:numel(ITEMS)*nSess)-1)/numel(ITEMS))+1);
                    catch E
                        if ~isempty(strfind(E.message,'sform and qform in the NIfTI header differ'))
                            ds=cosmo_fmri_dataset(fix_qform(fnTstats),'mask',fix_qform(fnMask),...
                                'targets',repmat(1:numel(ITEMS),1,nSess)','chunks',floor(((1:numel(ITEMS)*nSess)-1)/numel(ITEMS))+1);
                        else
                            aas_log(aap,true,E.message)
                        end
                    end
                    ds.sa.labels=cellfun(@(x) x{1}, repmat(ITEMS,1,nSess), 'UniformOutput', false)';
                    ds.sa.set=repmat((1:numel(ITEMS))',nSess,1);
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
            if cell_index(TASKS, 'RSA'), aas_log(aap,false,'INFO:Target DSM:'); disp(target_dsm); end
            
            %         imagesc(target_dsm)
            %         set(gca,'XTick',1:size(ds.samples,1),'XTickLabel',ds.sa.labels,...
            %             'YTick',1:size(ds.samples,1),'YTickLabel',ds.sa.labels)
            
            %% Run
            cosmo=cosmo_searchlight(ds,nbrhood,measure,measure_args);
            
            %         cosmo_plot_slices(ds_rsm_behav);
            
            % store results
            cosmo_fn=fullfile(TASKROOT,[t{1} 'map.nii']);
            cosmo_map2fmri(cosmo,cosmo_fn);
            
            %% Cleanup
            MVPA.unload;
            
            aap=aas_desc_outputs(aap,'subject',subj,[t{1} 'map'],cosmo_fn);
        end
    case 'checkrequirements'
        if ~aas_cache_get(aap,'cosmomvpa'), aas_log(aap,true,'CoSMoMVPA is not found'); end
        reqInps = {'mask' 'firstlevel_spm' 'firstlevel_spmts'};
        inps = aas_getstreams(aap,'input');
        missingInps = reqInps(cellfun(@(x) ~any(cell_index(inps,x)), reqInps));
        if ~isempty(missingInps), aas_log(aap,true,['Inputs not specified:' sprintf(' %s',missingInps{:})]); end
end
end

%% UTILS 
function val = check_qform(fnInput)
% Based on Chris Rorden's script in https://github.com/rordenlab/spmScripts at 25/11/2020

qOffsetBytes = 252;

fid = fopen(fnInput);
fseek(fid,qOffsetBytes,'bof');
qform_code = fread(fid,1,'int16');
fclose(fid);
val = qform_code ~= 0;
end

function fnOutput = fix_qform (fnInput)
% set qform to zero https://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1.h
% Based on Chris Rorden's script in https://github.com/rordenlab/spmScripts at 25/11/2020
qOffsetBytes = 252;
fnOutput = spm_file(fnInput,'prefix','z');

%create copy of image
copyfile(fnInput,fnOutput);

% modify qform
% - read file
fid = fopen(fnInput);
[data,count]=fread(fid, 'uint8');
fclose(fid);
% - modify both bytes of 16-bit qform_code
data(qOffsetBytes) = 0;
data(qOffsetBytes+1) = 0;
% - write file
fid = fopen(fnOutput,'w');
fwrite(fid,data);
fclose(fid);
end