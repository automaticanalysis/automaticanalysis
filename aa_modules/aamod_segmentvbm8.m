function [aap,resp] = aamod_segmentvbm8(aap, task, subjind)
% AAMOD_SEGMENT8 Perform VBM8's segment8 segmentation.
%
% [aap,resp] = AAMOD_SEGMENTVBM8(aap, task, subjind)
%
% These segmentations can then be fed into DARTEL or used in
% non-DARTEL VBM.
%
% input stream:     structural
% output streams:   seg8
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

resp='';

switch task
    case 'domain'
        resp='subject';
    case 'description'
        resp='VBM8 segment for structural images.'
    case 'doit'

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % options

        cfg.biasfwhm = aap.tasklist.currenttask.settings.biasfwhm; % 60
		cfg.biasreg = aap.tasklist.currenttask.settings.biasreg;   % .001
		cfg.affreg = aap.tasklist.currenttask.settings.affreg;     % 'mni';
        cfg.reg = aap.tasklist.currenttask.settings.reg;           % .001;
        cfg.vox = aap.tasklist.currenttask.settings.vox;           % voxel size things get resampled to
        cfg.mrf = aap.tasklist.currenttask.settings.mrf;           % markov random field cleanup
        cfg.samp = aap.tasklist.currenttask.settings.samp;         % sampling distance

        cfg.lkp = [1,1,2,2,3,3,4,4,4,5,5,5,5,6,6];
        cfg.writebiascorrected = aap.tasklist.currenttask.settings.writebiascorrected; % [bias_field bias_corrected]
        cfg.ngaus = aap.tasklist.currenttask.settings.ngaus;                           % [2 2 2 3 4 2];
        cfg.native = [1 1 0];                                                          % [native DARTEL_imported_rigid DARTEL_imported_affine]
        cfg.warped = [1 1 0];                                                          % normalised (don't write out for now, assume separate DARTEL)
        cfg.warpreg = 4;
        cfg.bb = {ones(2,3)*NaN};
        cfg.writedeffields = aap.tasklist.currenttask.settings.writedeffields;         % [1 1] would write them out

        cfg.cleanup = 0;
        cfg.dartelwarp = [];
        
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

        if ~strcmp(spm('ver'),'SPM8')
            aas_log(aap, 1, sprintf('%s requires SPM8.', mfilename));
        end

        if ~exist('spm_preproc_run', 'file')
            try
                % try adding a likely location
                addpath(fullfile(spm('dir'),'toolbox','Seg'))
            catch
            end
        end

        if ~exist('spm_preproc_run', 'file')
            aas_log(aap, true, 'spm_preproc8 is not in your Matlab path but needs to be.');
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


        % get the structural image
        img = aas_getfiles_bystream(aap, subjind, 'structural');


        if isempty(img) || strcmp(img,'/')
            aas_log(aap, true, 'Did not find a structural image.');
        end

        % if more than one found, use the first one and hope this is right
        img = strtok(img(1,:));

        aas_log(aap, false, sprintf('Found structural image: %s\n', img));


		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% do the segmentation

        tpm_nam = cfg.tpm;
        ngaus   = cfg.ngaus;
        nval    = {[1 0],[1 0],[1 0],[1 0],[1 0],[0 0]};
        for k=1:length(ngaus)
            tissue(k).tpm = [tpm_nam ',' num2str(k)]; % assign the tpm map
            tissue(k).ngaus = ngaus(k);  % and the number of gaussians
            tissue(k).native = cfg.native;
            tissue(k).warped = cfg.warped;
        end

        job.data = {img}; % VBM8 calls the image data rather than in job.channel
        
        job.channel(1).vols{1} = img;
        job.channel(1).biasreg = cfg.biasreg;
        job.channel(1).biasfwhm = cfg.biasfwhm;
        job.channel(1).write = cfg.writebiascorrected;
        job.channel(1).tpm = cfg.tpm;
        job.channel(1).ngaus = cfg.ngaus;
        job.channel(1).native = cfg.native;
        job.channel(1).warped = cfg.warped;

        job.tissue = tissue;

        job.warp.affreg = cfg.affreg;
        job.warp.reg = cfg.warpreg;
        job.warp.samp = cfg.samp;
        job.warp.write = cfg.writedeffields;
        job.warp.bb = cfg.bb;
        job.warp.vox = cfg.vox;
        job.warp.mrf = cfg.mrf;
        
        % these are for VBM8 options, defaults from cg_vbm8_defaults, but
        % these seem copied from the regular "new segment" so we use same
        % values from job.warp etc. above
        job.opts = [];
        job.opts.affreg = cfg.affreg;
        job.opts.samp = cfg.samp;
        job.opts.warpreg = cfg.warpreg;
        job.opts.tpm = {cfg.tpm};
        job.opts.ngaus = cfg.ngaus;
        job.opts.biasreg = cfg.biasreg;
        job.opts.biasfwhm = cfg.biasfwhm;
        
        
        job.extopts.sanlm = 1; % de-noising
        job.extopts.mrf = cfg.mrf;
        job.extopts.print = 0;
        job.extopts.cleanup = cfg.cleanup;
        job.extopts.dartelwarp = cfg.dartelwarp;
        
        job.output.warps = cfg.writedeffields;
        
        job.output.GM.native = cfg.native(1);
        job.output.GM.dartel=1;      % 1 = rigid aligned, 2 = affine aligned
        job.output.GM.warped = cfg.warped(1);
        job.output.GM.modulated = 1; % 1= normal SPM8-style modulation (but not relevant if we are not warping)
        
        job.output.WM.native = cfg.native(1);
        job.output.WM.dartel=1;      % 1 = rigid aligned, 2 = affine aligned
        job.output.WM.warped = cfg.warped(1);
        job.output.WM.modulated = 1; % 1= normal SPM8-style modulation (but not relevant if we are not warping)
        
        job.output.CSF.native = cfg.native(1);
        job.output.CSF.dartel=1;      % 1 = rigid aligned, 2 = affine aligned
        job.output.CSF.warped = cfg.warped(1);
        job.output.CSF.modulated = 1; % 1= normal SPM8-style modulation (but not relevant if we are not warping)
        
        job.output.bias.native = 0;
        job.output.bias.warped = 0;
        job.output.bias.affine = 0;
        
        job.output.label.native = 0;
        job.output.label.warped = 0;
        job.output.label.dartel = 0;
        
        job.output.jacobian.warped = 0;
        
        
        
        
        cg_vbm8_run(job);
        

        %% describe outputs

        [pth, nm, ext] = fileparts(img);
        seg8fn = fullfile(pth, sprintf('%s_seg8.mat', nm));
        aap = aas_desc_outputs(aap, subjind, 'seg8', seg8fn);

        tiss={'grey','white','csf'};
        for tissind=1:3
            aap = aas_desc_outputs(aap, subjind, sprintf('native_%s', tiss{tissind}), fullfile(pth, sprintf('p%d%s', tissind, [nm ext])));
            aap = aas_desc_outputs(aap, subjind, sprintf('dartelimported_%s', tiss{tissind}), fullfile(pth, sprintf('rp%d%s', tissind, [nm ext])));
            aap = aas_desc_outputs(aap, subjind, sprintf('normalised_density_%s', tiss{tissind}), fullfile(pth, sprintf('wp%d%s', tissind, [nm ext])));
            aap = aas_desc_outputs(aap, subjind, sprintf('normalised_volume_%s', tiss{tissind}), fullfile(pth, sprintf('mwp%d%s', tissind, [nm ext])));
        end
end