% AA module
% Runs BET (FSL Brain Extration Toolbox) on structural (usually)
% [For best functionality, it is recommended you run this after
% realignment and before writing the normalised EPI image
% If you do it before estimating the normalisation, make sure you normalise
% to a scull-stripped template, if at all possible!]

function [aap,resp]=aamod_fsl_FIRST(aap,task,subj)

resp='';

switch task
    case 'report'
        
    case 'doit'
        
        firstROIs =  {...
            10, 'L_Thal', 40; 49, 'R_Thal', 40; ...
            11, 'L_Caud', 30; 50, 'R_Caud', 30; ...
            12, 'L_Puta', 40; 51, 'R_Puta', 40; ...
            13, 'L_Pall', 40; 52, 'R_Pall', 40; ...
            16, 'BrStem', 40; ...
            17, 'L_Hipp', 30; 53, 'R_Hipp', 30; ...
            18, 'L_Amyg', 50; 54, 'R_Amyg', 50; ...
            26, 'L_Accu', 50; 58, 'R_Accu', 50};
        
        % Find out what stream we should BET
        inputstream = aap.tasklist.currenttask.inputstreams.stream;
        % And the names of the output streams
        outputstream = aap.tasklist.currenttask.outputstreams.stream;
        % And which are the streams which we output...
        outputstream = outputstream(~[strcmp(inputstream,outputstream)]);
        
        % Let us use the native space...
        Simg = aas_getfiles_bystream(aap,subj,inputstream{:});
        
        % Which file is considered, as determined by the structural parameter!
        if size(Simg,1) > 1
            Simg = deblank(Simg(aap.tasklist.currenttask.settings.structural, :));
            fprintf('WARNING: Several %s found, considering: \n', inputstream{:})
            for t = 1:length(aap.tasklist.currenttask.settings.structural)
                fprintf('\t%s\n', Simg(t,:))
            end
        end
        
        % Image that we will be using...
        [Spth Sfn Sext] = fileparts(Simg);
        
        % Normal, non-robust way of running FSL FIRST, will work for
        % low resolution, low contrast T1s
        if aap.tasklist.currenttask.settings.bet
            betOption = ' -b';
        else
            betOption = '';
        end
        
        if strcmp(aap.tasklist.currenttask.settings.mode, 'normal')
            %% Run FIRST
            fprintf('Running FSL FIRST\n')
            [junk, w]=aas_runfslcommand(aap, ...
                sprintf('run_first_all -v%s -i %s -o %s', ...
                betOption, ...
                Simg, Simg));
            
            V = spm_vol(fullfile(Spth, [Sfn '_all_fast_firstseg.nii']));
            Y = spm_read_vols(V);
            
            outSeg = {};
            
            %% Create masks for each structure
            for r = 1:length(firstROIs)
                M = Y==firstROIs{r,1};
                if sum(M(:)) > 0
                    V.fname = fullfile(Spth, [Sfn '_' firstROIs{r,2} Sext]);
                    outSeg = [outSeg, V.fname];
                    spm_write_vol(V,M);
                end
            end
        elseif strcmp(aap.tasklist.currenttask.settings.mode, 'weighted')
            %% 0) Create a mask to guide the FLIRT procedure
            mSimg = fullfile(Spth, ['m', Sfn, Sext]);
            copyfile(Simg, mSimg);
            img2mask(mSimg);            
            
            %% 1) FLIRT the subject brain to the template
            % FLIRT guide suggests to use as reference image the image with
            % larest contrast to noise and resolution
            % This is typically the structural and not the template...
            % FLIRTing Simg to Timg often goes wrong...
            
            FSLcommand = sprintf('first_flirt %s %s -inweight %s %s', ...
                Simg, ...
                fullfile(Spth, ['w' Sfn]), ...
                mSimg, ...
                betOption);
            aas_log(aap,false,FSLcommand);
            
            % Run command
            [junk, w] = aas_runfslcommand(aap, FSLcommand);
            
            % Check result
            spm_check_registration(strvcat(...
                Simg, ...
                fullfile(Spth, ['w', Sfn, Sext]), ...
                mSimg));
            
            %% 2) Run FIRST
            outSeg = {};
            for r = 1:size(firstROIs,1)
                ModelBMV = fullfile(aap.directory_conventions.fsldir, 'data', 'first', ...
                    'models_336_bin', [firstROIs{r, 2} '_bin.bmv']);
                if ~exist(ModelBMV, 'file')
                    ModelBMV = fullfile(aap.directory_conventions.fsldir, 'data', 'first', ...
                        'models_336_bin', '05mm', [firstROIs{r, 2} '_05mm.bmv']);
                end
                if ~exist(ModelBMV, 'file')
                    aas_log(aap, 1, 'Cannot find FSL model file!')
                end
                
                FSLcommand = sprintf('run_first -i %s -t %s -o %s -n %d -m %s', ...
                    Simg, ...
                    fullfile(Spth, ['w', Sfn, '.mat']), ...
                    fullfile(Spth, [Sfn '_' firstROIs{r, 2} Sext]), ...
                    firstROIs{r,3}, ...
                    ModelBMV);
                aas_log(aap,false,FSLcommand);
                
                outSeg = [outSeg, fullfile(Spth, [Sfn '_' firstROIs{r, 2} Sext])];
                
                % Run command
                [junk, w] = aas_runfslcommand(aap, FSLcommand);
            end
        elseif strcmp(aap.tasklist.currenttask.settings.mode, 'robust')
            %% 0) Get relevant template and mask
            
            % Find MNI template
            if aap.tasklist.currenttask.settings.bet
                Timg = deblank(ls(fullfile(aap.directory_conventions.fsldir, '*', 'standard', 'MNI152_T1_1mm_brain.*')));
            else
                Timg = deblank(ls(fullfile(aap.directory_conventions.fsldir, '*', 'standard', 'MNI152_T1_1mm.*')));
            end
            
            % Get the fileparts of the template
            [Tpth Tfn Text] = fileparts(Timg);
            if strcmp(Text, '.gz')
                [junk Tfn Text] = fileparts(Tfn);
                unix(['gunzip ' Timg])
                Timg = fullfile(Tpth, [Tfn Text]);
            end
            
            % Create a mask to guide the FLIRT procedure
            mSimg = fullfile(Spth, ['m', Sfn, Sext]);
            copyfile(Simg, mSimg);
            img2mask(mSimg);
                        
            % Find subcortical mask
            Mimg = deblank(ls(fullfile(aap.directory_conventions.fsldir, '*', 'standard', 'MNI152lin_T1_1mm_subbr_mask.*')));
            % Get the fileparts of the template
            [Mpth Mfn] = fileparts(Mimg);
            [junk Mfn Mext] = fileparts(Mfn);
            
            %% 1) FLIRT the template to our space
            % FLIRT guide suggests to use as reference image the image with
            % larest contrast to noise and resolution
            % This is typically the structural and not the template...
            % FLIRTing Simg to Timg often goes wrong...
            
            %FSLcommand = sprintf('flirt -in %s -ref %s -out %s -refweight %s -omat %s -cost %s', ...
            FSLcommand = sprintf('flirt -in %s -ref %s -out %s -omat %s -cost %s', ...
             Timg, ...
                Simg, ...
                fullfile(Spth, ['w', Tfn Text]), ... % template typically in nii.gz format...
                ... %mSimg, ...
                fullfile(Spth, [Tfn, '_1.mat']), ...
                aap.tasklist.currenttask.settings.cost);
            aas_log(aap,false,FSLcommand);
            
            % Run command
            [junk, w] = aas_runfslcommand(aap, FSLcommand);
            
            % Check result
            spm_check_registration(strvcat(...
                Simg, ...
                fullfile(Spth, ['w', Tfn Text])));
            
            %% 2) Set subcortical alignment mask into our space
            FSLcommand = sprintf('flirt -in %s -ref %s -out %s -applyxfm -init %s', ...
                Mimg, ...
                Simg, ...
                fullfile(Spth, ['w', Mfn Mext]), ... % mask typically in nii.gz format...
                fullfile(Spth, [Tfn, '_1.mat']));
            aas_log(aap,false,FSLcommand);
            
            % Run command
            [junk, w] = aas_runfslcommand(aap, FSLcommand);
            
            %% 3) Rerun FLIRT on warped template with subcortical mask...
            FSLcommand = sprintf('flirt -in %s -ref %s -out %s -omat %s -nosearch -refweight %s -cost %s', ...
                fullfile(Spth, ['w', Tfn Text]), ...
                Simg, ...
                fullfile(Spth, ['ww', Tfn Text]), ... % template typically in nii.gz format...
                fullfile(Spth, [Tfn, '_2.mat']), ...
                fullfile(Spth, ['w', Mfn Mext]), ...
                aap.tasklist.currenttask.settings.cost);
            aas_log(aap,false,FSLcommand);
            
            % Run command
            [junk, w] = aas_runfslcommand(aap, FSLcommand);
            
            % Check result
            spm_check_registration(strvcat(...
                Simg, ...
                fullfile(Spth, ['ww', Tfn Text])));
            
            %% 4) Invert .mat transforms and cocatenate them...
            FSLcommand = sprintf('convert_xfm -omat %s -inverse %s', ...
                fullfile(Spth, ['inv' Tfn, '_1.mat']), ...
                fullfile(Spth, [Tfn, '_1.mat']));
            aas_log(aap,false,FSLcommand);
            
            % Run command
            [junk, w] = aas_runfslcommand(aap, FSLcommand);
            
            FSLcommand = sprintf('convert_xfm -omat %s -inverse %s', ...
                fullfile(Spth, ['inv' Tfn, '_2.mat']), ...
                fullfile(Spth, [Tfn, '_2.mat']));
            aas_log(aap,false,FSLcommand);
            
            % Run command
            [junk, w] = aas_runfslcommand(aap, FSLcommand);
            
            FSLcommand = sprintf('convert_xfm -omat %s -concat %s %s', ...
                fullfile(Spth, ['inv' Tfn, '.mat']), ...
                fullfile(Spth, ['inv' Tfn, '_1.mat']), ...
                fullfile(Spth, ['inv' Tfn, '_2.mat']));
            aas_log(aap,false,FSLcommand);
            
            % Run command
            [junk, w] = aas_runfslcommand(aap, FSLcommand);
            
            %% 5) Run FIRST
            outSeg = {};
            for r = 1:size(firstROIs,1)
                ModelBMV = fullfile(aap.directory_conventions.fsldir, 'data', 'first', ...
                    'models_336_bin', [firstROIs{r, 2} '_bin.bmv']);
                if ~exist(ModelBMV, 'file')
                    ModelBMV = fullfile(aap.directory_conventions.fsldir, 'data', 'first', ...
                        'models_336_bin', '05mm', [firstROIs{r, 2} '_05mm.bmv']);
                end
                if ~exist(ModelBMV, 'file')
                    aas_log(aap, 1, 'Cannot find FSL model file!')
                end
                
                FSLcommand = sprintf('run_first -i %s -t %s -o %s -n %d -m %s', ...
                    Simg, ...
                    fullfile(Spth, ['inv' Tfn, '.mat']), ...
                    fullfile(Spth, [Sfn '_' firstROIs{r, 2} Sext]), ...
                    firstROIs{r,3}, ...
                    ModelBMV);
                aas_log(aap,false,FSLcommand);
                
                outSeg = [outSeg, fullfile(Spth, [Sfn '_' firstROIs{r, 2} Sext])];
                
                % Run command
                [junk, w] = aas_runfslcommand(aap, FSLcommand);
            end
        end
        
        %% DESCRIBE OUTPUTS!
        aap=aas_desc_outputs(aap,subj,'rois',outSeg);
        
        %% Save graphical output to common diagnostics directory
        subjname = aas_prepare_diagnostic(aap,subj);
        OVERcolours = aas_colours;
        
        %% Draw native template
        spm_check_registration(Simg)
        % Add segmentations...
        for t = 1:length(outSeg)
            spm_orthviews('addcolouredimage',1,outSeg{t}, OVERcolours{t})
        end
        
        spm_orthviews('reposition', [0 0 0])
        print('-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' subjname '.jpeg']));
        
        % Another diagnostic image, looking at how well the segmentation worked...
        Pthresh = 0.95;
        
        ROIdata = roi2hist(Simg, ...
            outSeg, Pthresh);
        
        [h, pv, ci, stats] = ttest2(ROIdata{2}, ROIdata{1});
        
        title(sprintf('GM vs WM... T-val: %0.2f (df = %d)', stats.tstat, stats.df))
        
        print('-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' subjname '_Hist.jpeg']));
        
        %% Diagnostic VIDEO
        if aap.tasklist.currenttask.settings.diagnostic
            Ydims = {'X', 'Y', 'Z'};
            
            for d = 1:length(Ydims)
                if (aap.tasklist.currenttask.settings.usesegmentnotnormalise)
                    aas_image_avi(Simg, ...
                        outSeg, ...
                        fullfile(aap.acq_details.root, 'diagnostics', [mfilename '__' subjname '_' Ydims{d} '.avi']), ...
                        d, ... % Axis
                        [800 600], ...
                        2, ... % Rotations
                        'none'); % No outline...
                    try close(2); catch; end
                end
            end
        end
end
