function [aap,resp] = aamod_segment8(aap, task, subjind)
% AAMOD_SEGMENT8 Perform SPM8's segment8 segmentation.
%
% [aap,resp] = AAMOD_SEGMENT8(aap, task, subjind)
%
% These segmentations can then be fed into DARTEL or used in
% non-DARTEL VBM.
%
% input stream:     structural
% output streams:   structural
%                    seg8
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
        resp='SPM8 segment8 for structural images.'
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
        cfg.native = [1 1];                                                            % [native DARTEL_imported]
        cfg.warped = [1 1];                                                            % normalised [modulated unmodulated]
        cfg.warpreg = 4;
        cfg.bb = {ones(2,3)*NaN};
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

        if job.warp.samp < 2
            fprintf('Note that the sampling distance is small, which means this might take quite a while (2-12+ hours depending on cluster load etc.)!\n');
        end

        spm_preproc_run(job);


        %% describe outputs
        aap = aas_desc_outputs(aap, subjind, 'structural', img);

        [pth, nm, ext] = fileparts(img);
        seg8fn = fullfile(pth, sprintf('%s_seg8.mat', nm));
        aap = aas_desc_outputs(aap, subjind, 'seg8', seg8fn);

        [pth nme ext]=fileparts(img);
        tiss={'grey','white','csf'};
        for tissind=1:3
            aap = aas_desc_outputs(aap, subjind, sprintf('native_%s', tiss{tissind}), fullfile(pth, sprintf('c%d%s', tissind, [nme ext])));
            aap = aas_desc_outputs(aap, subjind, sprintf('dartelimported_%s', tiss{tissind}), fullfile(pth, sprintf('rc%d%s', tissind, [nme ext])));
            aap = aas_desc_outputs(aap, subjind, sprintf('normalised_density_%s', tiss{tissind}), fullfile(pth, sprintf('wc%d%s', tissind, [nme ext])));
            aap = aas_desc_outputs(aap, subjind, sprintf('normalised_volume_%s', tiss{tissind}), fullfile(pth, sprintf('mwc%d%s', tissind, [nme ext])));
        end
end