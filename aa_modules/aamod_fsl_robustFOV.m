% AA module
% Runs BET (FSL Brain Extration Toolbox) on structural (usually)
% [For best functionality, it is recommended you run this after
% realignment and before writing the normalised EPI image
% If you do it before estimating the normalisation, make sure you normalise
% to a scull-stripped template, if at all possible!]

function [aap,resp]=aamod_fsl_robustFOV(aap,task,subj)

resp='';

switch task
    case 'report'
        
    case 'doit'
        
        % Find out what stream we should BET
        inputstream = aap.tasklist.currenttask.inputstreams.stream;
        % And the names of the output streams
        outputstream = aap.tasklist.currenttask.outputstreams.stream;
        % And which are the streams which we output...
        outputstream = outputstream(~[strcmp(inputstream,outputstream)]);
        
        % Get the image...
        Simg = aas_getfiles_bystream(aap,subj,inputstream{:});
        
        % Which file is considered, as determined by the structural parameter!        
        cSimg = deblank(Simg(aap.tasklist.currenttask.settings.structural, :));
        [Spth Sfn Sext] = fileparts(cSimg);
        if size(Simg,1) > 1
            fprintf('WARNING: Several %s found, considering: \n', inputstream{:})
            for t = 1:length(aap.tasklist.currenttask.settings.structural)
                fprintf('\t%s\n', Simg(t,:))
            end
        end        
        
        % Mat file of robust FOV...
        mSimg = fullfile(Spth, [Sfn '.mat']);
        imSimg = fullfile(Spth, ['inv' Sfn '.mat']);
        
        % Run robustFOV...
        fprintf('Running robustFOV\n')
        [junk, w]=aas_runfslcommand(aap, ...
            sprintf('robustfov -v -b %d -i %s -r %s -m %s', ...
            aap.tasklist.currenttask.settings.FOVsize, cSimg, cSimg, mSimg));
        
        % Get inverse transform?
        [junk, w]=aas_runfslcommand(aap, ...
            sprintf('convert_xfm -omat %s -inverse %s', imSimg, mSimg));        
        
        % Apply robustFOV on remaining images
        remImg = 1:size(Simg, 1);
        remImg = remImg(remImg ~= aap.tasklist.currenttask.settings.structural);
        for r = remImg
            [junk, w]=aas_runfslcommand(aap, ...
                sprintf('flirt -in %s -ref %s -applyxfm -init %s -out %s', ...
                deblank(Simg(r,:)), cSimg, imSimg, deblank(Simg(r,:))));
        end        
        
        %% DESCRIBE OUTPUTS!
        aap=aas_desc_outputs(aap,subj,inputstream{:},Simg);        
        
        %% DIAGNOSTIC IMAGE
        subjname = aas_prepare_diagnostic(aap,subj);
        
        %% Draw structural images...
        spm_check_registration(Simg)
        
        spm_orthviews('reposition', [0 0 0])
        
        print('-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' subjname '.jpeg']));
end
