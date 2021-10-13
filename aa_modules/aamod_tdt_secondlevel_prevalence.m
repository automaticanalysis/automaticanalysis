function [aap,resp]=aamod_tdt_secondlevel_prevalence(aap,task)
resp='';

switch task
    case 'report'
%         output = strrep(aas_getstreams(aap,'output'),'permuted_','');
%         fncfg = cellstr(spm_select('FPList',spm_file(aas_getfiles_bystream(aap,'subject',subj,aas_getstreams(aap,'output',1)),'path'),'.*cfg.mat$'));
%         dat = load(fncfg{1});
%         
%         if subj == 1 % first subject
%             for o = 1:numel(output)
%                 if  ~isfield(aap.report,sprintf('html_TDT%02d',o))
%                     aap.report.(sprintf('html_TDT%02d',o)).fname = fullfile(aap.report.condir,[aap.report.fbase sprintf('_TDT%02d.htm',o)]);
%                     aap = aas_report_add(aap,'C00',...
%                         sprintf('<a href="%s" target=_top>%s</a><br>',...
%                         aap.report.(sprintf('html_TDT%02d',o)).fname,...
%                         [strjoin(unique(dat.cfg.files.labelname,'stable'),' vs ') '; Measure: ' output{o}]));
%                     aap = aas_report_add(aap,sprintf('TDT%02d',o),['HEAD=' strjoin(unique(dat.cfg.files.labelname,'stable'),' vs ') '; Measure: ' output{o}]);
%                 end
%                 if ~isempty(aap.tasklist.currenttask.extraparameters.aap.directory_conventions.analysisid_suffix)
%                     aap = aas_report_add(aap,sprintf('TDT%02d',o),sprintf('<h2>Branch: %s</h2>',...
%                         aap.tasklist.currenttask.extraparameters.aap.directory_conventions.analysisid_suffix(2:end)));
%                 end
%             end
%         end
%          
%         for o = 1:numel(output)            
%             taskName = [strjoin(unique(dat.cfg.files.labelname,'stable'),' vs ') '; Measure: ' output{o}];
%             
%             aap = aas_report_add(aap,sprintf('TDT%02d',o),['Subject: ' basename(aas_getsubjpath(aap,subj)) '<br>']);
%             
%             aap = aas_report_add(aap,subj,sprintf('<h4>%02d. %s</h4>',o,taskName));
%             
%             imgfname = fullfile(aas_getsubjpath(aap,subj),...
%                 sprintf('diagnostic_aamod_tdt_firstlevel_statistics_%s_overlay_3_001.jpg',output{o}));
%             if exist(imgfname,'file')
%                 tstat = dlmread(strrep(imgfname,'_overlay_3_001.jpg','.txt'));
%                 
%                 % add image and report to subject
%                 aap = aas_report_add(aap, subj,'<table><tr>');
%                 aap = aas_report_add(aap, subj, sprintf('T = %2.2f - %2.2f</tr><tr>', tstat(1), tstat(2)));
%                 aap = aas_report_addimage(aap, subj, imgfname);
%                 aap = aas_report_add(aap,subj,'</tr></table>');
%                 
%                 % add image and report to summary
%                 aap = aas_report_add(aap,sprintf('TDT%02d',o),'<table><tr>');
%                 aap = aas_report_add(aap,sprintf('TDT%02d',o),sprintf('T = %2.2f - %2.2f</tr><tr>', tstat(1), tstat(2)));
%                 aap = aas_report_addimage(aap,sprintf('TDT%02d',o), imgfname);
%                 aap = aas_report_add(aap,sprintf('TDT%02d',o),'</tr></table>');
%             end
%         end
    case 'doit'
        [~, TDT] = aas_cache_get(aap,'tdt');
        TDT.load;
        
        statdir = fullfile(aas_getstudypath(aap),'stats');
        aas_makedir(aap,statdir);
        resfnameroot = fullfile(statdir, 'prevalence_');
        
        %% Background
        if ~isempty(aas_getsetting(aap,'overlay'))
            switch aas_getsetting(aap,'overlay.template')
                case 'averaged_structurals'
                    bgfname = [];
                    if aas_stream_has_contents(aap,'structural')
                        
                        for s = 1:numel(aap.acq_details.subjects)
                            bgfname{s} = aas_getfiles_bystream(aap,s,'structural');
                        end
                        aas_makedir(aap,fullfile(aas_getstudypath(aap),'structural'));
                        bgfname = spm_imcalc(char(bgfname),fullfile(aas_getstudypath(aap),'structural','background.nii'),'mean(X)',struct('dmtx',1,'mask',1));
                        bgfname = bgfname.fname;
                        
                    else
                        aas_log(aap, true, sprintf('%s: No structural stream found. Define one, or use SPMT1 template option. Exiting...', mfilename));
                    end
                    
                case 'SPMT1'
                    % assume a reasonable default location, but assume the user put
                    % the correct location in aap.directory_conventions.T1template if it's not empty
                    bgfname = 'toolbox/OldNorm/T1.nii';
                    if ~isempty(aap.directory_conventions.T1template), bgfname = aap.directory_conventions.T1template; end
                    if (bgfname(1) ~= '/'), bgfname = fullfile(fileparts(which('spm')),bgfname); end
                    if ~exist(bgfname,'file')
                        aas_log(aap, true, sprintf('%s: SPM T1 template not found. Exiting...', mfilename));
                    end
                    
                otherwise
                    aas_log(aap, true, sprintf('%s: Unknown template option. Exiting...', mfilename));
            end
        end
        
        dat = load(aas_getfiles_bystream(aap,'subject',1,'settings'));
        for decoding_measure = reshape(dat.cfg.results.output,1,[])
            %% Statistics
            resultfilenames = [resfnameroot decoding_measure{1}];
            
            % get data
            nSubj = aas_getN_bydomain(aap,'subject');
            inputimages = cell(nSubj,2);
            for subj = 1:nSubj
                res_image = cellstr(aas_getfiles_bystream(aap,'subject',subj,decoding_measure{1})); res_image = res_image(end);
                permuted_images = aas_getfiles_bystream(aap,'subject',subj,['permuted_' decoding_measure{1}]);
                % convert res_image into NIfTI if needed
                if strcmp(spm_file(res_image,'ext'),'mat')
                    V = spm_vol(permuted_images); V = V(1);
                    V.fname = spm_file(res_image{1},'ext','nii');
                    Y = read_image('matlab',read_header('matlab',res_image{1}));
                    spm_write_vol(V,Y);
                    res_image = cellstr(V.fname);
                end
                
                % get the original unpermuted result image as first image (required by the package)
                inputimages(subj, 1) = res_image;
                
                % put permuted images (4D) afterwards
                inputimages(subj, 2) = {permuted_images};
            end
            
            % run
            prevalenceTDT(inputimages, aas_getsetting(aap,'iteration'), resultfilenames);
            
            %% Overlay
            
            thrfnames = {};
            for f = cellstr(spm_select('FPList',statdir,['^' spm_file(resultfilenames,'basename') '_p[cu]{1}.*']))'
                % threshold
                V = spm_vol(f{1});
                Y = spm_read_vols(V);
                Y = 1-Y;
                Y = Y.*(Y>(1-aas_getsetting(aap,'threshold.p')));
                % % print global maximum
                % if any(Y,'all')
                %     [~,idx]=max(Y(:));
                %     indx=cell(ndims(Y),1);
                %     [indx{:}]=ind2sub(size(Y),idx); indx = vertcat(indx{:});
                %     fprintf('%s at\n',f{1});
                %     for l = 1:size(indx,2)
                %         fprintf('\t[%d %d %d]\n',indx(:,l))
                %     end
                % end
                V.fname = spm_file(f{1},'prefix','thresh_');
                spm_write_vol(V,Y);
                thrfnames{end+1} = V.fname;
                 
                if ~isempty(aas_getsetting(aap,'overlay'))
                    % edges of activation
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
                    % convert to mm
                    slims = sort(V.mat*slims,2);
                    % extend if too narrow (min. 50mm)
                    slims = slims + (repmat([-25 25],4,1).*repmat(diff(slims,[],2)<50,1,2));
                    
                    % - draw
                    axis = {'sagittal','coronal','axial'};
                    for a = 1:3
                        if any(Y,'all'), stat_fname = {V.fname}; else, stat_fname = {}; end
                        [fig, v] = map_overlay(bgfname,stat_fname,axis{a},slims(a,1):aas_getsetting(aap,'overlay.nth_slice'):slims(a,2));
                        fnsl{a} = fullfile(aas_getstudypath(aap), sprintf('diagnostic_aamod_tdt_firstlevel_statistics_%s_overlay_%d.jpg',spm_file(f{1},'basename'),a));
                        
                        if ~any(Y,'all')
                            annotation('textbox',[0 0.475 0.5 0.5],'String','No voxels survive threshold','FitBoxToText','on','fontweight','bold','color','y','fontsize',18,'backgroundcolor','k');
                        end
                        
                        spm_print(fnsl{a},fig,'jpg')
                        close(fig);
                    end
                    
                    dlmwrite(fullfile(aas_getstudypath(aap), sprintf('diagnostic_aamod_tdt_firstlevel_statistics_%s.txt',spm_file(f{1},'basename'))),[min(v(v~=0)), max(v)]);
                end
            end
            aap = aas_desc_outputs(aap,['thresholded_' decoding_measure{1}],char(thrfnames));
        end
        
        %% Cleanup
        TDT.unload;
        
    case 'checkrequirements'
        if ~aas_cache_get(aap,'tdt'), aas_log(aap,true,'TDT is not found'); end
        
        % automatic connection to the nearest aamod_tdt_decode
        srcmodulename = strrep(aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name,'secondlevel_prevalence','decode'); % CAVE: assumption on module name
        src = aas_getstreams(aas_setcurrenttask(aap,aas_getsourcestage(aap,srcmodulename,'settings')),'output'); src = setdiff(src,{'settings' 'mask'});
        inp = aas_getstreams(aap,'input'); 
        if any(strcmp(src,'settings_pairwise'))
            if ~any(strcmp(inp,'settings_pairwise')), aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'append','settings_pairwise','input'); end
            src = setdiff(src,{'settings_pairwise'});
        end
        inp = setdiff(inp,{'settings' 'settings_pairwise' 'structural'});
        for s = 1:numel(src) % every second excluding pairwise
            if strcmp(inp{s},src{s}), continue; end
            if s == 1
                aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'input',src{s},'input');
                aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'permuted_input',['permuted_' src{s}],'input');
                aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'thresholded_output',['thresholded_' src{s}],'output');                
            else
                aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'append',src{s},'input');
                aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'append',['permuted_' src{s}],'input');
                aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'append',['thresholded_' src{s}],'output');                
            end
            aas_log(aap,false,['INFO: ' aap.tasklist.currenttask.name ' input/output streams: ''' src{s} '''']);
        end
end
end