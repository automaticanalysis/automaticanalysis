function [aap,resp] = aamod_segment8(aap, task, subjind)
% AAMOD_SEGMENT8_MULTICHAN Perform SPM8's segment8 segmentation (multichannel).
%
% [aap,resp] = AAMOD_SEGMENT8(aap, task, subjind)
%
% These segmentations can then be fed into DARTEL or used in
% non-DARTEL VBM.
%
% input stream:     structural (may add others)
% output streams:   seg8
%                   forward_deformation
%                   inverse_deformation
%                   native_grey
%                   native_white
%                   native_csf
%                   dartelimported_grey
%                   dartelimported_white
%                   dartelimported_csf
%                   normalised_density_grey
%                   normalised_density_white
%                   normalised_density_csf
%                   normalised_volume_grey
%                   normalised_volume_white
%                   normalised_volume_csf
%
% A separate bias correction (e.g. with AAMOD_BIASCORRECT_SEGMENT8) prior to segmentation
% can improve the robustness.
%
% Multichannel functionality (jt 03/July/2012): Input images should be
% coregistered (e.g. with AAMOD_COREGISTER_GENERAL) prior to segmentation.
%
% This version also uses options writenormimg and startingaffineestimate
%  (jt 05/July/2012)
%

resp='';

switch task
    case 'report'
        fdiag = dir(fullfile(aas_getsubjpath(aap,subjind),'diagnostic_*.jpg'));
        if isempty(fdiag)
            diag(aap,subjind);
            fdiag = dir(fullfile(aas_getsubjpath(aap,subjind),'diagnostic_*.jpg'));
        end
        for d = 1:numel(fdiag)
            aap = aas_report_add(aap,subjind,'<table><tr><td>');
            imgpath = fullfile(aas_getsubjpath(aap,subjind),fdiag(d).name);
            aap=aas_report_addimage(aap,subjind,imgpath);
            [p f] = fileparts(imgpath); avipath = fullfile(p,[strrep(f(1:end-2),'slices','avi') '.avi']);
            if exist(avipath,'file'), aap=aas_report_addimage(aap,subjind,avipath); end
            aap = aas_report_add(aap,subjind,'</td></tr></table>');
        end
    case 'doit'

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % options

        cfg.biasfwhm = aap.tasklist.currenttask.settings.biasfwhm; % 60
		cfg.biasreg = aap.tasklist.currenttask.settings.biasreg;   % .001
		cfg.affreg = aap.tasklist.currenttask.settings.affreg;     % 'mni';
        cfg.reg = aap.tasklist.currenttask.settings.reg;           % .001;
        cfg.vox = aap.tasklist.currenttask.settings.vox;           % voxel size things get resampled to
        cfg.mrf = aap.tasklist.currenttask.settings.mrf;           % markov random field cleanup
        cfg.cleanup = aas_getsetting(aap,'cleanup'); if isempty(cfg.cleanup), cfg.cleanup = 0; end
        cfg.samp = aap.tasklist.currenttask.settings.samp;         % sampling distance

        cfg.lkp = [1,1,2,2,3,3,4,4,4,5,5,5,5,6,6];
        cfg.writebiascorrected = aap.tasklist.currenttask.settings.writebiascorrected; % [bias_field bias_corrected]
        cfg.ngaus = aap.tasklist.currenttask.settings.ngaus;                           % [2 2 2 3 4 2];
        cfg.native = [1 1];                                                            % [native DARTEL_imported]
        cfg.warped = [1 1];                                                            % normalised [modulated unmodulated]
        switch spm('ver')
            case 'SPM8'
                cfg.warpreg = 4;
                cfg.bb = {NaN(2,3)};                
            case {'SPM12b' 'SPM12'}
                cfg.warpreg = [0 1e-3 0.5 0.05 0.2];
                cfg.bb = NaN(2,3);
        end
        cfg.writedeffields = aap.tasklist.currenttask.settings.writedeffields;         % [1 1] would write them out

	    % If no full path to TPM specified, try to use the standard SPM one
        if isempty(aap.tasklist.currenttask.settings.tpm)
            cfg.tpm = fullfile(spm('dir'), 'toolbox', 'Seg', 'TPM.nii');
        else
            cfg.tpm = aap.tasklist.currenttask.settings.tpm;
        end

        if ~exist(cfg.tpm, 'file')
            aas_log(aap, true, sprintf('Specified TPM %s not found.', cfg.tpm));
        end

        aas_log(aap, false, sprintf('Segmenting using TPMs from %s.', cfg.tpm));


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % first make sure spm_preproc8 is in the path, as a toolbox it
        % might not be

        if ~exist('spm_preproc_run', 'file')
            switch spm('ver')
                case 'SPM8'
                    try
                        % try adding a likely location
                        addpath(fullfile(spm('dir'),'toolbox','Seg'))
                    catch
                    end
                case {'SPM12b' 'SPM12'}
                    try
                        % try adding a likely location
                        addpath(fullfile(spm('dir')));
                    catch
                    end                    
                otherwise
                    aas_log(aap, 1, sprintf('%s requires SPM8 or later.', mfilename));
            end
        end
        
        if ~exist('spm_preproc_run', 'file')
            aas_log(aap, true, 'spm_preproc_run is not in your Matlab path but needs to be.');
        end

        % make sure optimNn is in the path, usually with DARTEL
        if ~exist('optimNn', 'file')
            try
                addpath(fullfile(spm('dir'),'toolbox','DARTEL'))
            catch
            end
        end
        if ~exist('optimNn', 'file')
            aas_log(aap, true, 'optimNn is not in your Matlab path but needs to be.');
        end


        % get the images (multichan):
        channels = aap.tasklist.currenttask.settings.inputstreams.stream;
        if ~iscell(channels)
            channels={channels};
        end
		% check for empty streams
		ind = [];
		for c=1:length(channels)
            chaap = aap; chaap.options.verbose = -1;
            if aas_stream_has_contents(chaap,channels{c}) && ~isempty(aas_getfiles_bystream(chaap, subjind, channels{c}))
                ind = [ind c];
            end
		end
		channels = channels(ind);
		
        for c=1:length(channels)
            img{c} = aas_getfiles_bystream(aap, subjind, channels{c});

            if isempty(img{c}) || strcmp(img{c},'/')
                aas_log(aap, true, sprintf('Did not find a %s image.',channels{c}));
            end

            % if more than one found, use the first one and hope this is right
            img{c} = strtok(img{c}(1,:));
            
            aas_log(aap, false, sprintf('Found %s image: %s\n',channels{c},img{c}));
        end
        
        % For filenames:
        [pth, nm, ext] = fileparts(img{1});

        % Initial estimate of rigid-body rotation for 32-channel coil (jt 05/Jul/2012)
        StartingParameters=[0 0 0   0 0 0   1 1 1   0 0 0];
        ParameterFields={'x','y','z', 'pitch','roll','yaw', 'xscale','yscale','zscale', 'xaffign','yaffign','zaffign'};
        fnames=setdiff(fieldnames(aap.tasklist.currenttask.settings.affinestartingestimate),'COMMENT');
        for fieldind=1:length(fnames)
            % Which element in StartingParameters does this refer to?
        	whichitem=find([strcmp(fnames{fieldind},ParameterFields)]);
            % Generate a helpful error if it isn't recognised
            if (isempty(whichitem))
                err=sprintf('Unexpected field %s in header file aamod_segment8.xml - expected one of ',fnames{fieldind});
                err=[err sprintf('%s\t',ParameterFields)];
                aas_log(aap,true,err);
            end;
            % Put this in its place
            StartingParameters(whichitem)=aap.tasklist.currenttask.settings.affinestartingestimate.(fnames{fieldind});
        end;
            
        % Initial transform estimate to be passed to spm_preproc_run:
        StartingAffine = spm_matrix(StartingParameters);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% do the segmentation

        tpm_nam = cfg.tpm;
        ngaus   = cfg.ngaus;
        nval    = {[1 0],[1 0],[1 0],[1 0],[1 0],[0 0]};
        for k=1:length(ngaus)
            tissue(k).tpm    = [tpm_nam ',' num2str(k)]; % assign the tpm map
            tissue(k).ngaus  = ngaus(k);  % and the number of gaussians
            tissue(k).native = cfg.native;
            tissue(k).warped = cfg.warped;
        end
        
        % multichan (assumes same settings for all channels):
        for c=1:length(channels)
            job.channel(c).vols{1}  = img{c};
            job.channel(c).biasreg  = cfg.biasreg;
            job.channel(c).biasfwhm = cfg.biasfwhm;
            job.channel(c).write    = cfg.writebiascorrected;
            job.channel(c).tpm      = cfg.tpm;
            job.channel(c).ngaus    = cfg.ngaus;
            job.channel(c).native   = cfg.native;
            job.channel(c).warped   = cfg.warped;
        end

        job.tissue = tissue;

        job.warp.affreg         = cfg.affreg;
        job.warp.reg            = cfg.warpreg;
        job.warp.samp           = cfg.samp;
        job.warp.write          = cfg.writedeffields;
        job.warp.bb             = cfg.bb;
        job.warp.vox            = cfg.vox;
        job.warp.mrf            = cfg.mrf;
        job.warp.cleanup        = cfg.cleanup;
        job.warp.Affine         = StartingAffine; % jt
        
        if job.warp.samp < 2
            aas_log(aap,false,'Note that the sampling distance is small, which means this might take quite a while (2-12+ hours depending on cluster load etc.)!');
        end
        
        tic %(jt)
        aas_log(aap,false,sprintf('Starting to run segment8 job (%s)...', datestr(now, 'yyyy-mm-dd HH:MM')));
        spm_preproc_run(job);
        aas_log(aap,false,sprintf('\bDone in %.1f hours.', toc/60/60));
                

        % Find filenames ('c' imgs, y_ params named after diff channels?):
        for c=1:length(channels)
            [pth nm ext] = fileparts(img{c});
            if exist(fullfile(pth,['y_' nm ext]),'file');
                normparamfn=fullfile(pth,['y_' nm ext]);
                invparamfn=fullfile(pth,['iy_' nm ext]);
            end
            if exist(fullfile(pth,['c1' nm ext]),'file');
                c_img=img{c};
            end
        end
        
        % Apply deformation field to native structural? (jt 05/Jul/2012):
        if aap.tasklist.currenttask.settings.writenormimg
            opts = aap.tasklist.currenttask.settings.writenorm;
            for c=1:length(channels)
                aas_log(aap,false,'Applying normalisation parameters to input image(s)...');
                clear djob ojob
                ojob.ofname = '';
                ojob.fnames{1} = img{c};
                ojob.savedir.saveusr{1} = pth;
                ojob.interp = 1;
                switch spm('ver')
                    case 'SPM8'
                        spm_func_def = @spm_defs;
                        djob = ojob;
                        djob.comp{1}.def{1} = normparamfn;
                    case {'SPM12b' 'SPM12'}
                        spm_func_def = @spm_deformations;
                        ojob.mask = 0;
                        ojob.fwhm = opts.fwhm;                        
                        djob.comp{1}.def{1} = normparamfn;
                        djob.out{1}.savedef = ojob;
                        djob.out{2}.(opts.method) = ojob;
                        djob.out{2}.(opts.method).fov.bbvox.bb = cfg.bb;
                        djob.out{2}.(opts.method).fov.bbvox.vox = cfg.vox;
                        djob.out{2}.(opts.method).preserve = opts.preserve;
                    otherwise
                        aas_log(aap, 1, sprintf('%s requires SPM8 or later.', mfilename));
                end
                dout = spm_func_def(djob);
                aas_log(aap,false,'\b Done.');
            end
        end

        %% describe outputs

        [pth nm ext] = fileparts(c_img);
        seg8fn = fullfile(pth, sprintf('%s_seg8.mat', nm));
        aap = aas_desc_outputs(aap, subjind, 'seg8', seg8fn);

        tiss={'grey','white','csf'};
        for tissind=1:3
            aap = aas_desc_outputs(aap, subjind, sprintf('native_%s', tiss{tissind}), fullfile(pth, sprintf('c%d%s', tissind, [nm ext])));
            aap = aas_desc_outputs(aap, subjind, sprintf('dartelimported_%s', tiss{tissind}), fullfile(pth, sprintf('rc%d%s', tissind, [nm ext])));
            aap = aas_desc_outputs(aap, subjind, sprintf('normalised_density_%s', tiss{tissind}), fullfile(pth, sprintf('%sc%d%s', aap.spm.defaults.normalise.write.prefix, tissind, [nm ext])));
            aap = aas_desc_outputs(aap, subjind, sprintf('normalised_volume_%s', tiss{tissind}), fullfile(pth, sprintf('m%sc%d%s', aap.spm.defaults.normalise.write.prefix, tissind, [nm ext])));
        end
        
        % If user chose to write out normalised input image(s) (jt 05/Jul/2012)
        if aap.tasklist.currenttask.settings.writenormimg 
            pfx = aap.spm.defaults.normalise.write.prefix;
            if strcmp(opts.method,'push') && opts.preserve, pfx = ['m' pfx]; end
            if sum(opts.fwhm.^2)~=0, pfx = ['s' pfx]; end
            outstreams = aas_getstreams(aap,'output');
            for c=1:length(channels)
                aap = aas_desc_outputs(aap, subjind, outstreams{c}, fullfile(pth, sprintf('%s%s', pfx,[nm ext])));
            end
        end
        
        % If user chose to write out deformation fields (jt 05/Jul/2012)
        if aap.tasklist.currenttask.settings.writedeffields(1)
            aap = aas_desc_outputs(aap, subjind, 'forward_deformation_field', normparamfn);
        end
        if aap.tasklist.currenttask.settings.writedeffields(2)
            aap = aas_desc_outputs(aap, subjind, 'inverse_deformation_field', invparamfn);
        end
        
        %% Diagnostics
        if strcmp(aap.options.wheretoprocess,'localsingle')
            diag(aap,subjind);
        end
    case 'checkrequirements'
        %% Adjust outstream
        if ~aap.tasklist.currenttask.settings.writenormimg 
            for out = aap.tasklist.currenttask.settings.inputstreams.stream
                if any(strcmp(aas_getstreams(aap,'output'),out{1}))
                    aap = aas_renamestream(aap,aap.tasklist.currenttask.name,out{1},[],'output');
                    aas_log(aap,false,sprintf('REMOVED: %s output stream: %s', aap.tasklist.currenttask.name,out{1}));
                end
            end
        end
end
end

function diag(aap,subj) % [TA]
% SPM, AA
Simg = aas_getfiles_bystream(aap,subj,'structural');
localpath = aas_getsubjpath(aap,subj);
outSeg = char(...
    aas_getfiles_bystream(aap,subj,'native_csf'),...
    aas_getfiles_bystream(aap,subj,'native_white'),...
    aas_getfiles_bystream(aap,subj,'native_grey'));
outNSeg = char(...
    aas_getfiles_bystream(aap,subj,'normalised_density_csf'),...
    aas_getfiles_bystream(aap,subj,'normalised_density_white'),...
    aas_getfiles_bystream(aap,subj,'normalised_density_grey'));

OVERcolours = aas_colours;

%% Draw native template
spm_check_registration(Simg)
% Add normalised segmentations...
for r = 1:size(outSeg,1)
    spm_orthviews('addcolouredimage',1,deblank(outSeg(r,:)), OVERcolours{r});
end

spm_orthviews('reposition', [0 0 0])

try f = spm_figure('FindWin', 'Graphics'); catch; f = figure(1); end;
set(f,'Renderer','zbuffer');
print('-djpeg','-r150',...
    fullfile(localpath,['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_N_blob.jpg']));

%% Draw warped template
tmpfile = fullfile(aap.directory_conventions.fsldir,'data','standard','MNI152_T1_1mm.nii.gz'); % use FSL highres
gunzip(tmpfile,localpath); tmpfile = fullfile(localpath,'MNI152_T1_1mm.nii');

spm_check_registration(tmpfile)
% Add normalised segmentations...
for r = 1:size(outNSeg,1)
    spm_orthviews('addcolouredimage',1,outNSeg(r,:), OVERcolours{r})
end
spm_orthviews('reposition', [0 0 0])

try f = spm_figure('FindWin', 'Graphics'); catch; f = figure(1); end;
set(f,'Renderer','zbuffer');
print('-djpeg','-r150',...
    fullfile(localpath,['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_W_blob.jpg']));

close(f); clear global st;

%% Another diagnostic image, looking at how well the segmentation worked...
Pthresh = 0.95;

ROIdata = roi2hist(Simg, outSeg, Pthresh);

[h, pv, ci, stats] = ttest2(ROIdata{2}, ROIdata{1});

title(sprintf('GM vs WM... T(%d) = %0.2f, p = %1.4f', stats.df, stats.tstat, pv))

set(2,'Renderer','zbuffer');
print(2,'-djpeg','-r150',...
    fullfile(localpath,['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_Hist.jpg']));
try close(2); catch; end

%% Contours VIDEO of segmentations
aas_checkreg(aap,subj,outSeg(1,:),Simg)
aas_checkreg(aap,subj,outNSeg(1,:),tmpfile)

delete(tmpfile);
end