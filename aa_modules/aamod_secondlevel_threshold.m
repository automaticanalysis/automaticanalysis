% AA module - second level thresholding
% **********************************************************************
% You should no longer need to change this module - you may just
% modify the .xml or model in your user script
% **********************************************************************
% Tibor Auer MRC CBU Cambridge 2012-2013
%
% CHANGE HISTORY
%
% 04/2017 - fixed reporting. Now saves the results for each defined
% contrast instead of overwriting the same numel(SPM.xCon) files. Also
% replaced spaces with underscores in contrast names (because spaces
% in filenames is asking for trouble). [MSJ]
%
% 06/2018 - template selection is now explicit. Also removed 'continue'
% statements so even empty heatmaps are generated (missing results
% are confusing in the report). Added "no significant voxel"
% watermark to maps with no significant voxels. Add sanity check(s).
% Add optional "description" text which is overlayed on map if defined.
% [MSJ]
%

function [aap,resp]=aamod_secondlevel_threshold(aap,task)

resp='';

switch task
    
    case 'report'
        
        fnSPMs = aas_getfiles_bystream(aap, 'secondlevel_spm');
        
        for flc = 1:size(fnSPMs,1)
            
            fnSPM = deblank(fnSPMs(flc,:));
            loaded=load(fnSPM);
            SPM=loaded.SPM;
            
            [~,cname1,~] = fileparts(fileparts(fnSPM));
            
            for C = 1:numel(SPM.xCon)
                
                cname2 = SPM.xCon(C).name;
                
                aap = aas_report_add(aap,[],sprintf('<h4>%s %s</h4>', cname1, cname2));
                
                % recall we replaced spaces with underscores when we generated the diagnostic filename in 'doit'...
                
                cname2 = strrep(cname2,' ','_');
                
                % overlay_0 contains the axials -- we add that and the rendered
                % for each of the second level stats, which are probably
                % Group Mean (F), Mean Activation (>T), and Mean Deactivation (<T)
                
                f{1} = fullfile(aas_getstudypath(aap),...
                    sprintf('diagnostic_%s_%s_overlay_3_001.jpg', cname1, cname2));
                
                if exist(f{1},'file')
                    
                    tstat = dlmread(strrep(f{1},'_overlay_3_001.jpg','.txt'));
                    
                    f{2} = fullfile(aas_getstudypath(aap),...
                        sprintf('diagnostic_%s_%s_render.jpg', cname1, cname2));
                    
                    aap = aas_report_add(aap,[],'<table><tr>');
                    aap = aas_report_add(aap,[],sprintf('T = %2.2f - %2.2f</tr><tr>',tstat(1),tstat(2)));
                    
                    for i = 1:2
                        aap = aas_report_add(aap,[],'<td>');
                        aap=aas_report_addimage(aap,[],f{i});
                        aap = aas_report_add(aap,[],'</td>');
                    end
                    
                    aap = aas_report_add(aap,[],'</tr></table>');
                    
                end
                
            end
            
        end
        
        
    case 'doit'
        
        % sanity checks
        
        if (numel(aap.acq_details.subjects) < 4 && strcmp(aap.tasklist.currenttask.settings.threshold.correction,'FWE') )
            aas_log(aap, false, sprintf('WARNING: (%s): FWE may fail with fewer than 4 subjects.', mfilename));
        end
        
        if (strcmp(aap.tasklist.currenttask.settings.overlay.template, 'averaged_structurals'))
            aas_log(aap, false, sprintf('WARNING: (%s): structurals must be compatible (or normed) to use averaged_structurals template option', mfilename));
        end
        
        % Init
        u0   = aap.tasklist.currenttask.settings.threshold.p;				% height threshold
        nSl = aap.tasklist.currenttask.settings.overlay.nth_slice;
        tra = aap.tasklist.currenttask.settings.overlay.transparency;
        
        Outputs.thr = '';
        Outputs.cl = '';
        Outputs.sl = '';
        Outputs.Rend = '';
        
        cwd=pwd;
        localroot = aas_getstudypath(aap);
        
        nsub=length(aap.acq_details.subjects);
        aas_log(aap,false,sprintf('%d subjects',nsub));
        % New option to allow suffix to output file in extraparameters
        if (isfield(aap.tasklist.currenttask.extraparameters,'stats_suffix'))
            stats_suffix=aap.tasklist.currenttask.extraparameters.stats_suffix;
        else
            stats_suffix=[];
        end
        
        rfxrootdir = fullfile(aap.acq_details.root,[aap.directory_conventions.rfx stats_suffix]);
        cd(rfxrootdir);
        
        % template option is now explicit
        
        switch  aap.tasklist.currenttask.settings.overlay.template
            
            case 'averaged_structurals'
                
                tmpfile = [];
                
                inpstreams = aas_getstreams(aap,'input');
                
                if aas_stream_has_contents(aap,inpstreams{end})
                    
                    for s = 1:numel(aap.acq_details.subjects)
                        tmpfile{s} = aas_getfiles_bystream(aap,s,inpstreams{end});
                    end
                    aas_makedir(aap,fullfile(aas_getstudypath(aap),'structural'));
                    tmpfile = spm_imcalc(char(tmpfile),fullfile(aas_getstudypath(aap),'structural','background.nii'),'mean(X)',struct('dmtx',1,'mask',1));
                    tmpfile = tmpfile.fname;
                    
                else
                    aas_log(aap, true, sprintf('%s: No structural stream found. Define one, or use SPMT1 template option. Exiting...', mfilename));
                end
                
                
            case 'SPMT1'
                
                % assume a reasonable default location, but assume the user put
                % the correct location in aap.dir_con.T1template if it's not empty
                
                tmpfile = 'toolbox/OldNorm/T1.nii';
                if ~isempty(aap.directory_conventions.T1template) tmpfile = aap.directory_conventions.T1template; end
                if (tmpfile(1) ~= '/'), tmpfile = fullfile(fileparts(which('spm')),tmpfile); end
                
                if ~exist(tmpfile,'file')
                    aas_log(aap, true, sprintf('%s: SPM T1 template not found. Exiting...', mfilename));
                end
                
            otherwise
                
                aas_log(aap, true, sprintf('%s: Unknown template option. Exiting...', mfilename));
                
        end
        
        % Now get contrasts...
        
        fnSPMs = aas_getfiles_bystream(aap, 'secondlevel_spm');
        for flc = 1:size(fnSPMs,1)
            fnSPM = deblank(fnSPMs(flc,:));
            loaded=load(fnSPM);
            SPM=loaded.SPM;
            anadir=fileparts(fnSPM);
            cd(anadir);
            anadir = pwd; % resolve links
            fnSPM = spm_file(fnSPM,'path',anadir);
            
            for c = 1:numel(SPM.xCon)
                no_sig_voxels = false; % need this for later
                STAT = SPM.xCon(c).STAT;
                df = [SPM.xCon(c).eidf SPM.xX.erdf];
                XYZ  = SPM.xVol.XYZ;
                S    = SPM.xVol.S;   % Voxel
                R    = SPM.xVol.R;   % RESEL
                V = spm_vol(fullfile(anadir, SPM.xCon(c).Vspm.fname));
                Z = spm_get_data(SPM.xCon(c).Vspm,XYZ);
                dim = SPM.xCon(c).Vspm.dim;
                VspmSv   = cat(1,SPM.xCon(c).Vspm);
                n = 1; % No conjunction
                
                corr = aap.tasklist.currenttask.settings.threshold.correction;		% correction
                if strcmp(corr,'TFCE') && strcmp(STAT,'F'), corr = 'FWE'; end % TFCE does not support F-test -> FWE
                
                if strcmp(corr,'TFCE')
                    k = aas_getsetting(aap,'threshold.extent');
                    job.spmmat = {fnSPM};
                    job.mask = {fullfile(fileparts(fnSPM),'mask.nii,1')};
                    job.conspec = struct( ...
                        'titlestr','', ...
                        'contrasts',c, ...
                        'n_perm',5000, ...
                        'vFWHM',0 ...
                        );
                    job.openmp = 1;
                    cg_tfce_estimate(job);
                    iSPM = SPM;
                    iSPM.title = '';
                    iSPM.Ic = c;
                    iSPM.stattype = 'TFCE';
                    iSPM.thresDesc = corr;
                    iSPM.u = u0;
                    iSPM.k = k;
                    [SPM, xSPM] = cg_get_tfce_results(iSPM);
                    Z = xSPM.Z;
                    XYZ = xSPM.XYZ;
                    if isempty(Z)
                        aas_log(aap,false,sprintf('INFO: No voxels survive TFCE(%s)=%1.4f, k=%0.2g',corr, u0, k));
                        no_sig_voxels = true;
                    end
                else
                    % Height threshold filtering
                    switch corr
                        case 'iTT'
                            % TODO
                            [Z, XYZ, th] = spm_uc_iTT(Z,XYZ,u0,1);
                        case 'FWE'
                            u = spm_uc(u0,df,STAT,R,n,S);
                        case 'FDR'
                            u = spm_uc_FDR(u0,df,STAT,n,VspmSv,0);
                        case 'none'
                            u = spm_u(u0^(1/n),df,STAT);
                    end
                    Q      = find(Z > u);
                    Z      = Z(:,Q);
                    XYZ    = XYZ(:,Q);
                    if isempty(Q)
                        aas_log(aap,false,sprintf('INFO: No voxels survive height threshold u=%0.2g',u));
                        no_sig_voxels = true;
                    end
                    
                    % Extent threshold filtering
                    if ischar(aas_getsetting(aap,'threshold.extent')) % probability-based
                        k = strsplit(aas_getsetting(aap,'threshold.extent'),':'); k{2} = str2double(k{2});
                        iSPM = SPM;
                        iSPM.Ic = c;
                        iSPM.thresDesc = corr;
                        iSPM.u = u0;
                        iSPM.k = 0;
                        iSPM.Im = [];
                        [~,xSPM] = spm_getSPM(iSPM);
                        T = spm_list('Table',xSPM);
                        switch k{1}
                            case {'FWE' 'FDR'}
                                k{1} = ['p(' k{1} '-corr)'];
                            case {'none'}
                                k{1} = 'p(unc)';
                        end
                        pInd = strcmp(T.hdr(1,:),'cluster') & strcmp(T.hdr(2,:),k{1});
                        kInd = strcmp(T.hdr(2,:),'equivk');
                        k = min(cell2mat(T.dat(cellfun(@(p) ~isempty(p) && p<k{2}, T.dat(:,pInd)),kInd)));
                        if isempty(k), k = Inf; end
                    else
                        k = aas_getsetting(aap,'threshold.extent');
                    end
                    
                    A     = spm_clusters(XYZ);
                    Q     = [];
                    for i = 1:max(A)
                        j = find(A == i);
                        if length(j) >= k
                            Q = [Q j];
                        end
                    end
                    Z     = Z(:,Q);
                    XYZ   = XYZ(:,Q);
                    if isempty(Q)
                        aas_log(aap,false,sprintf('INFO: No voxels survive extent threshold k=%0.2g',k));
                        no_sig_voxels = true;
                    end
                end
                
                % Reconstruct
                Yepi  = zeros(dim(1),dim(2),dim(3));
                indx = sub2ind(dim,XYZ(1,:)',XYZ(2,:)',XYZ(3,:)');
                Yepi(indx) = Z;
                V.fname = spm_file(V.fname,'basename',strrep(spm_file(V.fname,'basename'),'spm','thr'));
                V.descrip = sprintf('thr{%s_%1.4f;ext_%d}%s',corr,u0,k,V.descrip(strfind(V.descrip,'}')+1:end));
                spm_write_vol(V,Yepi);
                
                % Cluster
                clusterfname = spm_file(V.fname,'prefix','cl_');
                if ~isempty(aas_getsetting(aap,'cluster'))
                    if no_sig_voxels
                        copyfile(V.fname,clusterfname);
                    else
                        switch aas_getsetting(aap,'cluster.method')
                            case 'fusionwatershed'
                                [s,FWS] = aas_cache_get(aap,'fws');
                                if ~s
                                    aas_log(aap,false,'WARNING: Fusion-Watershed is not installed! --> clustering skipped');
                                else
                                    FWS.load;
                                    settings = aas_getsetting(aap,'cluster.options.fusionwatershed');
                                    obj = fws.generate_ROI(V.fname,...
                                        'threshold_method','z','threshold_value',0.1,...
                                        'filter',settings.extentprethreshold,'radius',settings.searchradius,'merge',settings.mergethreshold,...
                                        'plot',false,'output',true);
                                    
                                    % - exclude small (<k) ROIs
                                    smallROIs = obj.table.ROIid(obj.table.Volume < k);
                                    obj.label(reshape(arrayfun(@(l) any(l == smallROIs), obj.label(:)), obj.grid.d)) = 0;
                                    obj.table(arrayfun(@(l) any(l == smallROIs), obj.table.ROIid),:) = [];
                                    
                                    % - save results
                                    save.vol(obj.label,obj.grid,spm_file(clusterfname,'ext',''),'Compressed',false);
                                    writetable(obj.table,spm_file(clusterfname,'ext','csv'));
                                    FWS.unload;
                                end
                        end
                    end
                end
                
                % Overlay
                % - edges of activation (in mm)
                slims = ones(4,2);
                sAct = arrayfun(@(x) any(Yepi(x,:,:),'all'), 1:size(Yepi,1));
                if numel(find(sAct))<2, slims(1,:) = [1 size(Yepi,1)];
                else, slims(1,:) = [find(sAct,1,'first') find(sAct,1,'last')]; end
                sAct = arrayfun(@(y) any(Yepi(:,y,:),'all'), 1:size(Yepi,2));
                if numel(find(sAct))<2, slims(2,:) = [1 size(Yepi,2)];
                else, slims(2,:) = [find(sAct,1,'first') find(sAct,1,'last')]; end
                sAct = arrayfun(@(z) any(Yepi(:,:,z),'all'), 1:size(Yepi,3));
                if numel(find(sAct))<2, slims(3,:) = [1 size(Yepi,3)];
                else, slims(3,:) = [find(sAct,1,'first') find(sAct,1,'last')]; end
                % - convert to mm
                slims = sort(V.mat*slims,2);
                % - extend if too narrow (min. 50mm)
                slims = slims + (repmat([-25 25],4,1).*repmat(diff(slims,[],2)<50,1,2));
                
                % - draw
                axis = {'sagittal','coronal','axial'};
                for a = 1:3
                    if ~no_sig_voxels, stat_fname = {V.fname}; else, stat_fname = {}; end
                    [fig, v] = map_overlay(tmpfile,stat_fname,axis{a},slims(a,1):nSl:slims(a,2));
                    if (~isempty(aap.tasklist.currenttask.settings.description))
                        annotation('textbox',[0 0.5 0.5 0.5],'String',aap.tasklist.currenttask.settings.description,'FitBoxToText','on','fontweight','bold','color','y','fontsize',18,'backgroundcolor','k');
                    end
                    if (no_sig_voxels)
                        annotation('textbox',[0 0.475 0.5 0.5],'String','No voxels survive threshold','FitBoxToText','on','fontweight','bold','color','y','fontsize',18,'backgroundcolor','k');
                    end
                    
                    % strip the first level contrast name from the SPM.mat directory
                    [~,cname1,~] = fileparts(fileparts(fnSPM));
                    
                    % strip the second level contrast from the SPM structure
                    % this will probably be 'Group Mean' (i.e., an F test), 'Group Mean
                    % Activation (t-test) or 'Group Mean Deactivation' (t-test). We swap
                    % any spaces for underscores bc spaces in filesnames is asking for
                    % Unix trouble [MSJ]
                    cname2 = SPM.xCon(c).name; cname2 = strrep(cname2, ' ', '_');
                    fnsl{a} = fullfile(localroot, sprintf('diagnostic_%s_%s_overlay_%d.jpg', cname1, cname2, a));
                    
                    if (~isempty(aap.tasklist.currenttask.settings.description))
                        annotation('textbox',[0 0.5 0.5 0.5],'String',aap.tasklist.currenttask.settings.description,'FitBoxToText','on','fontweight','bold','color','y','fontsize',18,'backgroundcolor','k');
                    end
                    
                    if (no_sig_voxels)
                        annotation('textbox',[0 0.475 0.5 0.5],'String','No voxels survive threshold','FitBoxToText','on','fontweight','bold','color','y','fontsize',18,'backgroundcolor','k');
                    end
                    
                    fnsl{a} = spm_print(fnsl{a},fig,'jpg');
                end
                
                dlmwrite(fullfile(localroot, sprintf('diagnostic_%s_%s.txt', cname1, cname2)),[min(v(v~=0)), max(v)]);
                
                % Render
                
                % this fails if no sig voxels
                
                if ~no_sig_voxels
                    if numel(Z)  < 2
                        Z = horzcat(Z,Z);
                        XYZ = horzcat(XYZ,XYZ);
                    end
                end
                
                dat.XYZ = XYZ;
                dat.t = Z';
                dat.mat = SPM.xVol.M;
                dat.dim = dim;
                rendfile  = aap.directory_conventions.Render;
                if ~exist(rendfile,'file') && (rendfile(1) ~= '/'), rendfile = fullfile(fileparts(which('spm')),rendfile); end
                fn3d = fullfile(localroot, sprintf('diagnostic_%s_%s_render.jpg', cname1, cname2));
                global prevrend
                prevrend = struct('rendfile',rendfile, 'brt',0.5, 'col',eye(3));
                out = spm_render(dat,0.5,rendfile); spm_figure('Close','Graphics');
                img = vertcat(horzcat(out{1},out{3},out{5}),horzcat(out{2},out{4},out{6}));
                fig = figure;
                imshow(img,'Border','tight');
                
                if (~isempty(aap.tasklist.currenttask.settings.description))
                    annotation('textbox',[0 0.5 0.5 0.5],'String',aap.tasklist.currenttask.settings.description,'FitBoxToText','on','fontweight','bold','color','y','fontsize',18,'backgroundcolor','k');
                end
                
                if (no_sig_voxels)
                    annotation('textbox',[0 0.45 0.5 0.5],'String','No voxels survive threshold','FitBoxToText','on','fontweight','bold','color','y','fontsize',18,'backgroundcolor','k');
                end
                
                print(fig,'-noui',fn3d,'-djpeg','-r300');
                close(fig);
                
                % Outputs
                if exist(V.fname,'file'), Outputs.thr = strvcat(Outputs.thr, V.fname); end
                if exist(clusterfname,'file'), Outputs.cl = strvcat(Outputs.cl, clusterfname); end                
                for f = 1:numel(fnsl)
                    if exist(fnsl{f},'file'), Outputs.sl = strvcat(Outputs.sl, fnsl{f}); end
                end
                if exist(fn3d,'file'), Outputs.Rend = strvcat(Outputs.Rend, fn3d); end
                clear fnsl
            end
        end
        
        cd (cwd);
        
        % Describe outputs
        
        aap=aas_desc_outputs(aap,'secondlevel_thr',Outputs.thr);
        aap=aas_desc_outputs(aap,'secondlevel_clusters',Outputs.cl);
        aap=aas_desc_outputs(aap,'secondlevel_thrslice',Outputs.sl);
        aap=aas_desc_outputs(aap,'secondlevel_thr3D',Outputs.Rend);        
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end

end