% AA module - first level thresholding
% **********************************************************************
% You should no longer need to change this module - you may just
% modify the .xml or model in your user script
% **********************************************************************
% Tibor Auer MRC CBU Cambridge 2012-2013

function [aap,resp]=aamod_firstlevel_threshold(aap,task,subj)

resp='';

switch task
    case 'report'
        % Collect all contrast names and prepare summary
        cons = [aap.tasksettings.aamod_firstlevel_contrasts(aap.tasklist.currenttask.index).contrasts(2:end).con];
        conNames = {cons.name};
        [junk, a] = unique(conNames,'first');
        conNames = conNames(sort(a));
        if subj == 1 % first
            for C = 1:numel(conNames)
                if  ~isfield(aap.report,sprintf('html_C%02d',C))
                    aap.report.(sprintf('html_C%02d',C)).fname = fullfile(aap.report.condir,[aap.report.fbase sprintf('_C%02d.htm',C)]);
                    aap = aas_report_add(aap,'C00',...
                        sprintf('<a href="%s" target=_top>%s</a><br>',...
                        aap.report.(sprintf('html_C%02d',C)).fname,...
                        ['Contrast: ' conNames{C}]));
                    aap = aas_report_add(aap,sprintf('C%02d',C),['HEAD=Contrast: ' conNames{C}]);
                end
                if ~isempty(aap.tasklist.currenttask.extraparameters.aap.directory_conventions.analysisid_suffix)
                    aap = aas_report_add(aap,sprintf('C%02d',C),sprintf('<h2>Branch: %s</h2>',...
                        aap.tasklist.currenttask.extraparameters.aap.directory_conventions.analysisid_suffix(2:end)));
                end
            end
        end
        
        % Now get contrasts...
        fSPM = aas_getfiles_bystream(aap, subj,'firstlevel_spm');
        loaded = load(fSPM);
        
        for C = 1:numel(loaded.SPM.xCon)
            conName = strrep_multi(loaded.SPM.xCon(C).name,{' ' ':' '>'},{'' '_' '-'});
            conInd = find(strcmp(conNames, loaded.SPM.xCon(C).name));
            
            if isempty(conInd), continue; end % implicit contrasts (e.g. eachagainstbaseline)
            
            % Study summary
            aap = aas_report_add(aap,sprintf('C%02d',conInd),['Subject: ' basename(aas_getsubjpath(aap,subj)) '<br>']);
            
            % Single subject
            aap = aas_report_add(aap,subj,sprintf('<h4>%02d. %s</h4>',conInd,conName));
            f{1} = fullfile(aas_getsubjpath(aap,subj),...
                sprintf('diagnostic_aamod_firstlevel_threshold_C%02d_%s_overlay_0.jpg',conInd,conName));
            if exist(f{1},'file')
                tstat = dlmread(strrep(f{1},'_overlay_0.jpg','.txt'));
                f{2} = fullfile(aas_getsubjpath(aap,subj),...
                    sprintf('diagnostic_aamod_firstlevel_threshold_C%02d_%s_render.jpg',conInd,conName));
                
                % Add images to Single subject report
                aap = aas_report_add(aap,subj,'<table><tr>');
                aap = aas_report_add(aap,subj,sprintf('T = %2.2f - %2.2f</tr><tr>',tstat(1),tstat(2)));
                for i = 1:2
                    aap = aas_report_add(aap,subj,'<td>');
                    aap=aas_report_addimage(aap,subj,f{i});
                    aap = aas_report_add(aap,subj,'</td>');
                end
                aap = aas_report_add(aap,subj,'</tr></table>');
                
                % Add images to Study summary
                aap = aas_report_add(aap,sprintf('C%02d',conInd),'<table><tr>');
                aap = aas_report_add(aap,sprintf('C%02d',conInd),sprintf('T = %2.2f - %2.2f</tr><tr>',tstat(1),tstat(2)));
                for i = 1:2
                    aap = aas_report_add(aap,sprintf('C%02d',conInd),'<td>');
                    aap=aas_report_addimage(aap,sprintf('C%02d',conInd),f{i});
                    aap = aas_report_add(aap,sprintf('C%02d',conInd),'</td>');
                end
                aap = aas_report_add(aap,sprintf('C%02d',conInd),'</tr></table>');
                                
            end
        end
        
    case 'doit'
        % Init
        try doTFCE = aap.tasklist.currenttask.settings.threshold.doTFCE; catch, doTFCE = 0; end % TFCE?
        corr = aap.tasklist.currenttask.settings.threshold.correction;  % correction
        u0   = aap.tasklist.currenttask.settings.threshold.p;            % height threshold
        k   = aap.tasklist.currenttask.settings.threshold.extent;       % extent threshold {voxels}
        nSl = aap.tasklist.currenttask.settings.overlay.nth_slice;
        tra = aap.tasklist.currenttask.settings.overlay.transparency;
        Outputs.thr = '';
        Outputs.sl = '';
        Outputs.Rend = '';
        
        cwd=pwd;
        localroot = aas_getsubjpath(aap,subj);
        
        % get the subdirectories in the main directory
        anadir = fullfile(localroot,aap.directory_conventions.stats_singlesubj);
        cd(anadir);
        
        if aas_stream_has_contents(aap,subj,'structural') % Structural if available (backward compatibility)
            tmpfile = aas_getfiles_bystream(aap, subj,'structural');
            if size(tmpfile,1) > 1 % in case of norm_write (first: native, second: normalised)
                tmpfile = tmpfile(2,:); 
            end
        else  % Template
            aas_log(aap,false,'Structural cannot be loaded! Template will be used...');
            tmpfile = aap.directory_conventions.T1template;
            if tmpfile(1) ~= '/'
                if exist(tmpfile,'file'), tmpfile = which(tmpfile);
                else tmpfile = fullfile(fileparts(which('spm')),tmpfile); end
            end
        end
                
        Vtemplate=spm_vol(tmpfile);
        [Ytemplate, tXYZ]=spm_read_vols(Vtemplate);
        tXYZ=[tXYZ;ones(1,size(tXYZ,2))];
        % work out threshold for template
        threshprop=0.10;
        ys=sort(Ytemplate(~isnan(Ytemplate)));
        bright3=ys(round(length(ys)*0.3));
        bright97=ys(round(length(ys)*0.97));
        thresh=bright3*(1-threshprop)+bright97*threshprop;
        Ytemplate=Ytemplate.*(Ytemplate>thresh);
        
        % Now get contrasts...
        SPM=[]; 
        fSPM = aas_getfiles_bystream(aap, subj,'firstlevel_spm');
        load(fSPM);
        
        for c = 1:numel(SPM.xCon)
            conName = strrep_multi(SPM.xCon(c).name,{' ' ':' '>'},{'' '_' '-'});
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

            if doTFCE
                job.spmmat = {fSPM};
                job.mask = {fullfile(fileparts(fSPM),'mask.nii,1')};
                job.conspec = struct( ...
                    'titlestr','', ...
                    'contrasts',c, ...
                    'n_perm',5000, ...
                    'vFWHM',0 ...
                    );
                job.tbss = 0;
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
                    aas_log(aap,false,sprintf('WARNING: No voxels survive TFCE(%s)=%1.4f, k=%0.2g',corr, u0, k));
                    continue;
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
                    aas_log(aap,false,sprintf('WARNING: No voxels survive height threshold u=%0.2g',u));
                    continue;
                end
                
                % Extent threshold filtering
                A     = spm_clusters(XYZ);
                Q     = [];
                for i = 1:max(A)
                    j = find(A == i);
                    if length(j) >= k;
                        Q = [Q j];
                    end
                end
                Z     = Z(:,Q);
                XYZ   = XYZ(:,Q);
                if isempty(Q)
                    aas_log(aap,false,sprintf('WARNING: No voxels survive extent threshold k=%0.2g',k));
                    continue;
                end
            end
            
            % Reconstruct
            Yepi  = zeros(dim(1),dim(2),dim(3));
            indx = sub2ind(dim,XYZ(1,:)',XYZ(2,:)',XYZ(3,:)');
            Yepi(indx) = Z;
            vname = spm_file(V.fname,'basename');
            V.fname = spm_file(V.fname,'basename',strrep(vname,'spm','thr'));
            V.descrip = sprintf('thr{%s_%1.4f;ext_%d}%s',corr,u0,k,V.descrip(strfind(V.descrip,'}')+1:end));
            spm_write_vol(V,Yepi);
            
            % Resize
            iXYZ=SPM.xVol.M\tXYZ;
            rYepi=spm_sample_vol(Yepi,iXYZ(1,:),iXYZ(2,:),iXYZ(3,:),0);
            rYepi=reshape(rYepi,size(Ytemplate));
            
            % Overlay
            fnsl = '';
            for a = 0:2 % in 3 axes
                arYepi = shiftdim(rYepi,a);
                aYtemplate = shiftdim(Ytemplate,a);
                
                for iSl = 1:nSl % Adjust slice selection according to the activation
                    iYepi = arYepi(:,:,iSl:nSl:end);
                    if any(iYepi(:)~=0), break; end
                end
                iYepi = img_rot90(iYepi);
                iYtemplate = img_rot90(aYtemplate(:,:,iSl:nSl:end));
                
                [img, cm, v] = map_overlay(iYtemplate,iYepi,1-tra);                                
                mon = tr_3Dto2D(img_tr(img(:,:,:,1),a==2));
                mon(:,:,2) = tr_3Dto2D(img_tr(img(:,:,:,2),a==2));
                mon(:,:,3) = tr_3Dto2D(img_tr(img(:,:,:,3),a==2));
                fnsl(a+1,:) = fullfile(localroot, sprintf('diagnostic_aamod_firstlevel_threshold_C%02d_%s_overlay_%d.jpg',c,conName,a));
                imwrite(mon,deblank(fnsl(a+1,:)));
            end
            dlmwrite(strrep(deblank(fnsl(1,:)),'_overlay_0.jpg','.txt'),[min(v(v~=0)), max(v)]);
            % Render
            if numel(Z)  < 2 % render fails with only one active voxel
                Z = horzcat(Z,Z);
                XYZ = horzcat(XYZ,XYZ);
            end
            % render fails with single first slice 
            for a = 1:3
                if all(XYZ(a,:)==1)
                    Z = horzcat(Z,Z(end));
                    XYZ = horzcat(XYZ,XYZ(:,end)+circshift([1;0;0],a-1));
                end
            end
            dat.XYZ = XYZ;
            dat.t = Z';
            dat.mat = SPM.xVol.M;
            dat.dim = dim;
            rendfile  = aap.directory_conventions.Render;
            if ~exist(rendfile,'file') && (rendfile(1) ~= '/'), rendfile = fullfile(fileparts(which('spm')),rendfile); end
            fn3d = fullfile(localroot,sprintf('diagnostic_aamod_firstlevel_threshold_C%02d_%s_render.jpg',c,conName));
            global prevrend
            prevrend = struct('rendfile',rendfile, 'brt',0.5, 'col',eye(3));
            out = spm_render(dat,0.5,rendfile); spm_figure('Close','Graphics');
            clear img; for i = 1:numel(out), img(1:size(out{i},1),1:size(out{i},2),:,i) = out{i}; end
            mon = tr_3Dto2D(squeeze(img(:,:,1,[1 3 5 2 4 6])));
            mon(:,:,2) = tr_3Dto2D(squeeze(img(:,:,2,[1 3 5 2 4 6])));
            mon(:,:,3) = tr_3Dto2D(squeeze(img(:,:,3,[1 3 5 2 4 6])));
            mon = mon(1:size(mon,2)*2/3,:,:);
            imwrite(mon,fn3d);
            
            % Outputs
            if exist(V.fname,'file'), Outputs.thr = strvcat(Outputs.thr, V.fname); end
            for f = 1:size(fnsl,1)
                if exist(fnsl(f,:),'file'), Outputs.sl = strvcat(Outputs.sl, fnsl(f,:)); end
            end
            if exist(fn3d,'file'), Outputs.Rend = strvcat(Outputs.Rend, fn3d); end
        end
        cd (cwd);
        
        % Describe outputs
        aap=aas_desc_outputs(aap,subj,'firstlevel_thr',Outputs.thr);
        aap=aas_desc_outputs(aap,subj,'firstlevel_thrslice',Outputs.sl);
        aap=aas_desc_outputs(aap,subj,'firstlevel_thr3D',Outputs.Rend);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
end

%% UTILS

function fo = img_rot90(fi)
for i = 1:size(fi,3)
    fo(:,:,i) = rot90(fi(:,:,i),1);
end
end

function fo = img_tr(fi,toDo)
if nargin < 2, toDo = true; end
if toDo
    nslice = size(fi,3);
    for i = 1:nslice
        fo(:,:,i) = fliplr(rot90(fi(:,:,i),1));
    end
else
    fo = fi;
end
end