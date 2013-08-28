function [aap, resp] = aamod_biascorrect_segment8(aap, task, subjind)
%AAMOD_BIASCORRECT_SEGMENT8 Bias corrects a structural image using SPM8's "new segment".
%
% [aap, resp] = AAMOD_BIASCORRECT_SEGMENT8(aap, task, subjind) Performs the bias
% correction for the first structural iamge found for this subject.
%
% Performing a separate bias correction can be helpful for tissue class segmentation,
% particularly with multi-channel head coils, which can introduce significant
% inhomogeneities. A typical processing pipeline would be to bias correct the raw T1 image
% and then pass the bias-corrected T1 image to segmentation.
%
% input stream:     structural
%
% output streams:   structural  the bias-corrected t1 image
%                   seg8        the seg8.mat file generated from segmentation

resp = '';

switch task
	case 'domain'
		resp = 'subject';
	case 'description'
		resp = 'SPM8 segment8 to bias correct structural image';
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
        cfg.writebiascorrected = [0 1];                            % [bias_field bias_corrected]
        cfg.ngaus = aap.tasklist.currenttask.settings.ngaus;       % [2 2 2 3 4 2];
        cfg.native = [0 0];                                        % native and DARTEL imported
        cfg.warped = [0 0];                                        % normalised modulated (not unmodulated)
        cfg.warpreg = 4;
        cfg.bb = {ones(2,3)*NaN};
        cfg.writedeffields = [0 0];                                % [1 1] would write them out


        % If no full path to tissue probability map (TPM) specified, try to use the standard SPM one
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
        %img = aas_getfiles_bystream(aap, subjind, 'structural');
        inStream = aap.tasklist.currenttask.inputstreams.stream{1};
        img = aas_getfiles_bystream(aap, subjind, inStream);
        

        if isempty(img) || strcmp(img,'/')
            aas_log(aap, true, 'Did not find a structural image.');
        end

        % if more than one found, use the first one and hope this is right
        img = strtok(img(1,:));

        aas_log(aap, false, sprintf('Found structural image: %s\n', img));


		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% do the segmentation (and bias correction)

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

        % get the filename for the bias-corrected image (which has 'm' prepended)
        [pth, nm, ext] = fileparts(img);
        mimgfn = fullfile(pth, sprintf('m%s%s', nm, ext));

        % the bias-corrected structural image replaces the input image in the stream
        aap = aas_desc_outputs(aap, subjind, inStream, mimgfn);

        % get file name for *seg8.mat file
        seg8fn = fullfile(pth, sprintf('%s_seg8.mat', nm));
        aap = aas_desc_outputs(aap, subjind, 'seg8', seg8fn);
end