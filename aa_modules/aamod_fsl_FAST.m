% AA module
% Runs BET (FSL Brain Extration Toolbox) on structural (usually)
% [For best functionality, it is recommended you run this after
% realignment and before writing the normalised EPI image
% If you do it before estimating the normalisation, make sure you normalise
% to a scull-stripped template, if at all possible!]

function [aap,resp]=aamod_fsl_FAST(aap,task,subj)

resp='';

switch task
    case 'report'
        
    case 'doit'
        
        segMapping = [3 1 2];
        
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
            aas_log(aap,false,sprintf('WARNING: Several %s found, considering: ', inputstream{:}))
            for t = 1:length(aap.tasklist.currenttask.settings.structural)
                aas_log(aap,false,sprintf('\t%s', Simg(t,:)))
            end
        end
        
        % Image that we will be using for BET...
        [Spth Sfn Sext] = fileparts(Simg);
        
        % Run BET [-R Using robust setting to improve performance!]
        aas_log(aap,false,'Running FSL FAST')
        [junk, w]=aas_runfslcommand(aap, ...
            sprintf('fast %s %s ', ...
            aap.tasklist.currenttask.settings.options, Simg));
               
        %% FIND OUTPUT
        % Get the segmented images
        D = dir(fullfile(Spth, '*pve_*'));
        % Move these images to something that looks like the SPM output
        outSeg = cell(1,length(segMapping));
        
        for d = 1:length(segMapping)
            outSeg{d} = fullfile(Spth, ['c' num2str(d) Sfn Sext]);
            unix(['mv ' fullfile(Spth, D(segMapping(d)).name) ' ' outSeg{d}]);
        end
        
        %% DESCRIBE OUTPUTS!
        aap=aas_desc_outputs(aap,subj,'segmentation',outSeg);
        
        %% Save graphical output to common diagnostics directory
        subjname = aas_prepare_diagnostic(aap,subj);
        
        % This will only work for 1-7 segmentations
        OVERcolours = aas_colours;
        
        %% Draw native template
        spm_check_registration(Simg)
        % Add segmentations...
        for t = 1:length(segMapping)
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
