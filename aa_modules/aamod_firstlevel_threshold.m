% AA module - first level thresholding
% **********************************************************************
% You should no longer need to change this module - you may just
% modify the .xml or model in your user script
% **********************************************************************
% Based on spm_getSPM, aas_overlaywithtemplate
% Tibor Auer MRC CBU Cambridge 2012-2013

function [aap,resp]=aamod_firstlevel_threshold(aap,task,subj)

resp='';

switch task
    case 'domain'
        resp='subject';   % this module needs to be run once per subject
        
    case 'report'
        CON = aap.tasksettings.aamod_firstlevel_contrasts.contrasts(2).con;
        for C = 1:numel(CON)
            % Study summary
            if ~isfield(aap.report,sprintf('html_C%02d',C)) % new contrast
                aap.report.(sprintf('html_C%02d',C)).fname = fullfile(aap.report.condir,[aap.report.fbase sprintf('_C%02d.htm',C)]);
                aap = aas_report_add(aap,'C00',...
                    sprintf('<a href="%s" target=_top>%s</a><br>',...
                    aap.report.(sprintf('html_C%02d',C)).fname,...
                    ['Contrast: ' CON(C).name]));
                aap = aas_report_add(aap,sprintf('C%02d',C),['HEAD=Contrast: ' CON(C).name]);                
            end
            aap = aas_report_add(aap,sprintf('C%02d',C),['Subject: ' basename(aas_getsubjpath(aap,subj)) '<br>']);
            
            % Single subject
            aap = aas_report_add(aap,subj,sprintf('<h4>%02d. %s</h4>',C,CON(C).name));
            f{1} = fullfile(aas_getsubjpath(aap,subj),...
                sprintf('diagnostic_aamod_firstlevel_threshold_C%02d_%s_overlay.jpg',C,CON(C).name));
            if exist(f{1},'file')
                f{2} = fullfile(aas_getsubjpath(aap,subj),...
                    sprintf('diagnostic_aamod_firstlevel_threshold_C%02d_%s_render.jpg',C,CON(C).name));
                
                % Add images to Single subject report
                aap = aas_report_add(aap,subj,'<table><tr>');
                for i = 1:2
                    aap = aas_report_add(aap,subj,'<td>');
                    aap=aas_report_addimage(aap,subj,f{i});
                    aap = aas_report_add(aap,subj,'</td>');
                end
                aap = aas_report_add(aap,subj,'</tr></table>');
                
                % Add images to Study summary
                aap = aas_report_add(aap,sprintf('C%02d',C),'<table><tr>');
                for i = 1:2
                    aap = aas_report_add(aap,sprintf('C%02d',C),'<td>');
                    aap=aas_report_addimage(aap,sprintf('C%02d',C),f{i});
                    aap = aas_report_add(aap,sprintf('C%02d',C),'</td>');
                end
                aap = aas_report_add(aap,sprintf('C%02d',C),'</tr></table>');
                                
            end
        end
        
    case 'doit'
        
        cwd=pwd;
        % get the subdirectories in the main directory
        subj_dir = aas_getsubjpath(aap,subj);
        
        anadir = fullfile(subj_dir,aap.directory_conventions.stats_singlesubj);
        cd(anadir);
        
        % Init
        corr = aap.tasklist.currenttask.settings.threshold.correction;  % correction
        u0   = aap.tasklist.currenttask.settings.threshold.p;            % height threshold
        k   = aap.tasklist.currenttask.settings.threshold.extent;       % extent threshold {voxels}
        nSl = aap.tasklist.currenttask.settings.overlay.nth_slice;
        tra = aap.tasklist.currenttask.settings.overlay.transparency;
        Outputs.thr = '';
        Outputs.sl = '';
        Outputs.Rend = '';
        
        % Template
        tmpfile = aap.directory_conventions.T1template;
        if ~exist(tmpfile,'file') && (tmpfile(1) ~= '/'), tmpfile = fullfile(fileparts(which('spm')),tmpfile); end
        Vtemplate=spm_vol(tmpfile);
        [Ytemplate tXYZ]=spm_read_vols(Vtemplate);
        tXYZ=[tXYZ;ones(1,size(tXYZ,2))];
        % work out threshold for template
        threshprop=0.10;
        ys=sort(Ytemplate(~isnan(Ytemplate)));
        bright3=ys(round(length(ys)*0.3));
        bright97=ys(round(length(ys)*0.97));
        thresh=bright3*(1-threshprop)+bright97*threshprop;
        Ytemplate=Ytemplate.*(Ytemplate>thresh);
        
        % Collect first-level spms
        cType = {aap.tasksettings.aamod_firstlevel_contrasts.contrasts(2).con.type};
        con_dir = aas_getsubjpath(aap,subj,cell_index({aap.tasklist.main.module.name},'aamod_firstlevel_contrast'));
        con_dir=fullfile(con_dir,aap.directory_conventions.stats_singlesubj);
        if cell_index(cType,'T')
            copyfile(fullfile(con_dir,'spmT*'),anadir);
        end
        if cell_index(cType,'F')
            copyfile(fullfile(con_dir,'spmF*'),anadir);
        end
        
        % Now get contrasts...
        SPM=load(aas_getfiles_bystream(aap,subj,'firstlevel_spm')); SPM = SPM.SPM;
        for c = 1:numel(SPM.xCon)
            STAT = SPM.xCon(c).STAT;
            df = [SPM.xCon(c).eidf SPM.xX.erdf];
            XYZ  = SPM.xVol.XYZ;
            S    = SPM.xVol.S;   % Voxel
            R    = SPM.xVol.R;   % RESEL
            V = spm_vol(SPM.xCon(c).Vspm.fname);
            Z = spm_get_data(SPM.xCon(c).Vspm,XYZ);
            dim = SPM.xCon(c).Vspm.dim;
            VspmSv   = cat(1,SPM.xCon(c).Vspm);
            n = 1; % No conjunction
            
            % Height threshold filtering
            switch corr
                case 'iTT'
                    % TODO
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
                fprintf('\n');
                warning('No voxels survive height threshold u=%0.2g',u);
                continue;
            end
            
            % Extent threshold filtering
            A     = spm_clusters(XYZ);
            Q     = [];
            for i = 1:max(A)
                j = find(A == i);
                if length(j) >= k; Q = [Q j]; end
            end
            Z     = Z(:,Q);
            XYZ   = XYZ(:,Q);
            if isempty(Q)
                fprintf('\n');
                warning('No voxels survive extent threshold k=%0.2g',k);
                continue;
            end
            
            % Reconstruct
            Yepi  = zeros(dim(1),dim(2),dim(3));
            indx = sub2ind(dim,XYZ(1,:)',XYZ(2,:)',XYZ(3,:)');
            Yepi(indx) = Z;
            V.fname = strrep(V.fname,'spm','thr');
            V.descrip = sprintf('thr{%s_%1.4f;ext_%d}%s',corr,u0,k,V.descrip(findstr(V.descrip,'}')+1:end));
            spm_write_vol(V,Yepi);
            
            % Resize
            iXYZ=SPM.xVol.M\tXYZ;
            rYepi=spm_sample_vol(Yepi,iXYZ(1,:),iXYZ(2,:),iXYZ(3,:),0);
            rYepi=reshape(rYepi,size(Ytemplate));
            
            % Overlay
            for iSl = 1:nSl % Adjust slice selection according to the activation
                iYepi = rYepi(:,:,iSl:nSl:end);
                if any(iYepi(:)~=0), break; end
            end
            iYepi = img_rot90(iYepi);
            iYtemplate = img_rot90(Ytemplate(:,:,iSl:nSl:end));
            [img cm v] = map_overlay(iYtemplate,iYepi,1-tra);
            mon = tr_RGBtoMontage(img);
            f = figure;
            montage(mon);
            colormap(cm);
            cb = colorbar;
            yT = [walley(v) find(v<0, 1, 'last' ) find(v>0, 1 ) peak(v)];
            set(cb,'YTick',yT,'YTickLabel',v(yT));
            fnsl = fullfile(aas_getsubjpath(aap,subj), sprintf('diagnostic_aamod_firstlevel_threshold_C%02d_%s_overlay.jpg',c,SPM.xCon(c).name));
            print(f,'-djpeg','-r150',fnsl);
            
            close(f);
            
            % Render
            if numel(Z)  < 2 % Render fails with only one active voxel
                Z = horzcat(Z,Z);
                XYZ = horzcat(XYZ,XYZ);
            end
            dat.XYZ = XYZ;
            dat.t = Z';
            dat.mat = SPM.xVol.M;
            dat.dim = dim;
            rendfile  = aap.directory_conventions.Render;
            if ~exist(rendfile,'file') && (rendfile(1) ~= '/'), rendfile = fullfile(fileparts(which('spm')),rendfile); end
            spm_render(dat,0.5,rendfile);
            fn3d = fullfile(aas_getsubjpath(aap,subj),sprintf('diagnostic_aamod_firstlevel_threshold_C%02d_%s_render.jpg',c,SPM.xCon(c).name));
            saveas(spm_figure('GetWin','Graphics'),fn3d);
            img = imread(fn3d); img = img(size(img,1)/2:end,:,:); imwrite(img,fn3d);
            
            % Outputs
            if exist(fullfile(anadir,V.fname),'file'), Outputs.thr = strvcat(Outputs.thr,fullfile(anadir,V.fname)); end
            if exist(fnsl,'file'), Outputs.sl = strvcat(Outputs.sl,fnsl); end
            if exist(fn3d,'file'), Outputs.Rend = strvcat(Outputs.Rend,fn3d); end
            
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