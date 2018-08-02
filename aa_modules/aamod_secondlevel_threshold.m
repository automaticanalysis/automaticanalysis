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

function [aap,resp]=aamod_secondlevel_threshold(aap,task,subj)

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
 					sprintf('diagnostic_%s_%s_overlay_0.jpg', cname1, cname2));
		
                if exist(f{1},'file')
					
                    tstat = dlmread(strrep(f{1},'_overlay_0.jpg','.txt'));
					
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

        % Init
		
        try doTFCE = aap.tasklist.currenttask.settings.threshold.doTFCE; catch, doTFCE = 0; end % TFCE?
        corr = aap.tasklist.currenttask.settings.threshold.correction;		% correction
        u0   = aap.tasklist.currenttask.settings.threshold.p;				% height threshold
        k   = aap.tasklist.currenttask.settings.threshold.extent;			% extent threshold {voxels}
        nSl = aap.tasklist.currenttask.settings.overlay.nth_slice;
        tra = aap.tasklist.currenttask.settings.overlay.transparency;
		
        Outputs.thr = '';
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
        end;
        
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
						if size(tmpfile{s},1) > 1 % in case of norm_write (first: native, second: normalised)
							tmpfile{s} = tmpfile{s}(2,:);
						end
					end

				end
				
				% if we didn't find a structural stream, try to load the SPM template rather than just halting 
				% (this to be consistent with how previous versions of the module worked...)
				
				if isempty(tmpfile)
					tmpfile = 'toolbox/OldNorm/T1.nii';
					if ~isempty(aap.directory_conventions.T1template) tmpfile = aap.directory_conventions.T1template; end
					if (tmpfile(1) ~= '/'), tmpfile = fullfile(fileparts(which('spm')),tmpfile); end	
					if exist(tmpfile,'file')
						aas_log(aap, false, sprintf('WARNING (%s): Structural not found. Using SPM T1 template instead...', mfilename));
					else
						aas_log(aap, true, sprintf('%s: Cannot find normalised_structural. Exiting...', mfilename));
					end
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
		
        
        Vtemplate=spm_vol(tmpfile); if iscell(Vtemplate), Vtemplate = cell2mat(Vtemplate); end
        [Ytemplate, tXYZ]=spm_read_vols(Vtemplate); if ndims(Ytemplate) == 4, Ytemplate = mean(Ytemplate,4); end
        tXYZ=[tXYZ;ones(1,size(tXYZ,2))];
        % work out threshold for template
        threshprop=0.10;
        ys=sort(Ytemplate(~isnan(Ytemplate)));
        bright3=ys(round(length(ys)*0.3));
        bright97=ys(round(length(ys)*0.97));
        thresh=bright3*(1-threshprop)+bright97*threshprop;
        Ytemplate=Ytemplate.*(Ytemplate>thresh);
        
        % Now get contrasts...
        fnSPMs = aas_getfiles_bystream(aap, 'secondlevel_spm');
        for flc = 1:size(fnSPMs,1)
            fnSPM = deblank(fnSPMs(flc,:));
            loaded=load(fnSPM);
            SPM=loaded.SPM;
            anadir=fileparts(fnSPM);
            cd(anadir);
            
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
                
                if doTFCE
                    job.spmmat = {fSPM};
                    job.mask = {fullfile(fileparts(fSPM),'mask.nii,1')};
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
						aas_log(aap,false,sprintf('INFO: No voxels survive extent threshold k=%0.2g',k));
						no_sig_voxels = true;
					end
                end
                
                % Reconstruct
                Yepi  = zeros(dim(1),dim(2),dim(3));
                indx = sub2ind(dim,XYZ(1,:)',XYZ(2,:)',XYZ(3,:)');
                Yepi(indx) = Z;
                V.fname = strrep(V.fname,'spm','thr');
                V.descrip = sprintf('thr{%s_%1.4f;ext_%d}%s',corr,u0,k,V.descrip(strfind(V.descrip,'}')+1:end));
                spm_write_vol(V,Yepi);
                
                % Resize
                iXYZ=SPM.xVol.M\tXYZ;
                rYepi=spm_sample_vol(Yepi,iXYZ(1,:),iXYZ(2,:),iXYZ(3,:),0);
                rYepi=reshape(rYepi,size(Ytemplate));
                
                % Overlay
				
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
					if (~isempty(aap.tasklist.currenttask.settings.description))
						mon = insertInImage(mon, @()text(40,25,aap.tasklist.currenttask.settings.description),...
						{'fontweight','bold','color','y','fontsize',18,...
						'linewidth',1,'margin',5,'backgroundcolor','k'});	
					end
					if (no_sig_voxels)
						mon = insertInImage(mon, @()text(40,50,'No voxels survive threshold'),...
						{'fontweight','bold','color','y','fontsize',18,...
						'linewidth',1,'margin',5,'backgroundcolor','k'});	
					end	

					% strip the first level contrast name from the SPM.mat directory

					[~,cname1,~] = fileparts(fileparts(fnSPM));

					% strip the second level contrast from the SPM structure
					% this will probably be 'Group Mean' (i.e., an F test), 'Group Mean
					% Activation (t-test) or 'Group Mean Deactivation' (t-test). We swap
					% any spaces for underscores bc spaces in filesnames is asking for
					% Unix trouble [MSJ]

					cname2 = SPM.xCon(c).name; cname2 = strrep(cname2, ' ', '_');
					fnsl(a+1,:) = fullfile(localroot, sprintf('diagnostic_%s_%s_overlay_%d.jpg', cname1, cname2, a));
					imwrite(mon,deblank(fnsl(a+1,:)));
					
				end
				
                dlmwrite(strrep(deblank(fnsl(1,:)),'_overlay_0.jpg','.txt'),[min(v(v~=0)), max(v)]);
				
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
                clear img; for i = 1:numel(out), img(1:size(out{i},1),1:size(out{i},2),:,i) = out{i}; end
                mon = tr_3Dto2D(squeeze(img(:,:,1,[1 3 5 2 4 6])));
                mon(:,:,2) = tr_3Dto2D(squeeze(img(:,:,2,[1 3 5 2 4 6])));
                mon(:,:,3) = tr_3Dto2D(squeeze(img(:,:,3,[1 3 5 2 4 6])));
                mon = mon(1:size(mon,2)*2/3,:,:);
				if (~isempty(aap.tasklist.currenttask.settings.description))
					mon = insertInImage(mon, @()text(40,25,aap.tasklist.currenttask.settings.description),...
					{'fontweight','bold','color','y','fontsize',18,...
					'linewidth',1,'margin',5,'backgroundcolor','k'});	
				end
				if (no_sig_voxels)
					mon = insertInImage(mon, @()text(40,50,'No voxels survive threshold'),...
					{'fontweight','bold','color','y','fontsize',18,...
					'linewidth',1,'margin',5,'backgroundcolor','k'});	
				end	
                imwrite(mon,fn3d);
                
                % Outputs
                if exist(V.fname,'file'), Outputs.thr = strvcat(Outputs.thr, V.fname); end
                for f = 1:size(fnsl,1)
                    if exist(fnsl(f,:),'file'), Outputs.sl = strvcat(Outputs.sl, fnsl(f,:)); end
                end
                if exist(fn3d,'file'), Outputs.Rend = strvcat(Outputs.Rend, fn3d); end
                clear fnsl
            end
        end
        
        cd (cwd);
        
        % Describe outputs
		
        aap=aas_desc_outputs(aap,'secondlevel_thr',Outputs.thr);
        aap=aas_desc_outputs(aap,'secondlevel_thrslice',Outputs.sl);
        aap=aas_desc_outputs(aap,'secondlevel_thr3D',Outputs.Rend);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end

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

