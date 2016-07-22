% AA module
% Use the niftyseg (& niftyreg) toolboxes to segment the brain
% ...and optionally make a brain Mask (MATLAB Image Manipulation Toolbox required)

function [aap,resp]=aamod_LoAd(aap,task,subj)

resp='';

switch task
    case 'doit'
        
        % Set LoAd paths...
        pth=getenv('PATH');
        setenv('PATH',[pth ...
            ':' fullfile(aap.directory_conventions.niftysegdir, 'bin') ...
            ':' fullfile(aap.directory_conventions.niftyregdir, 'bin') ...
            ]);
        
        % Set libraries
        libpth = getenv('LD_LIBRARY_PATH');
        setenv('LD_LIBRARY_PATH', [libpth ...
            ':' fullfile(aap.directory_conventions.niftysegdir, 'lib') ...
            ':' fullfile(aap.directory_conventions.niftyregdir, 'lib') ...
            ]);
        
        % Set libraries (mac)
        libMpth = getenv('DYLD_LIBRARY_PATH');
        setenv('DYLD_LIBRARY_PATH', [libMpth ...
            ':' fullfile(aap.directory_conventions.niftysegdir, 'lib') ...
            ':' fullfile(aap.directory_conventions.niftyregdir, 'lib') ...
            ]);
        
        %% Get structural & mask
        
        Simg = aas_getfiles_bystream(aap,subj,'structural');
        % Which file is considered, as determined by the structural parameter!
        if size(Simg,1) > 1
            Simg = deblank(Simg(aap.tasklist.currenttask.settings.structural, :));
            aas_log(aap,false,sprintf('WARNING: Several structurals found, considering:\n\t%s', Simg(1,:)))
        end
        [Spth, Sfn, Sext] = fileparts(Simg);
        
        BETmask=aas_getfiles_bystream(aap,subj,'BETmask');
        for b = 1:size(BETmask,1)
            if ~isempty(strfind(deblank(BETmask(b,:)), 'brain_mask.nii'))
                BETmask = deblank(BETmask(b,:));
                break
            elseif b == size(BETmask,1)
                aas_log(aap, true, 'Cannot find the correct mask')
            end
        end
        
        % Let us loop the load algorhythm several times, so as to obtain
        % a single WM and GM mass...
        
        aas_log(aap,false,'Runing LoAd...')
        %% Use LoAd to segment the structural!
        LoAd_command = ['sh LoAd_brainonly.sh ' ... % Run LoAd command
            Simg ' ' ... % structural
            BETmask]; % mask
        
        [s w] = aas_shell(LoAd_command);
        
        %% Use seg_maths to extract the relevant
        outSeg = '';
        
        tissues = {'WM' 'GM' 'CSF' 'DeepGM' 'Ventricles'};
        
        for t = 1:length(tissues)
            % For each tissue...
            Mfn = fullfile(Spth, [tissues{t} Sext]);
            outSeg = strvcat(outSeg, Mfn);
            segmaths_command = ['seg_maths ' ... % Segment the end result...
                fullfile(Spth, [Sfn '_segmentation' Sext]) ' ' ... % Segmentation image
                '-tp ' num2str(t-1) ' ' ...
                Mfn];
            [s w] = aas_shell(segmaths_command);
        end
        
        %% BET mask
        mY = 0;
        % There are 5 sensible tissue classes, the rest are not
        % Exclude CSF, as this will make the brain mask too large...
        for t = [1 2 4 5]
            mY = mY + spm_read_vols(spm_vol(deblank(outSeg(t,:))));
        end
        mY = mY > 0;
        
        
        % Fill any holes that may be remaining
        try
            mY = imfill(mY,'holes');
        catch aa_error
            aas_log(aap, false, 'Could not fill the holes in brain mask')
        end
        
        V = spm_vol(BETmask);
        BETmask = fullfile(Spth, [Sfn '_LoADbrain_mask' Sext ]);
        V.fname = BETmask;
        spm_write_vol(V,double(mY));
        
        spm_check_registration(strvcat( ...
            Simg, ... % Get structural
            BETmask)); % Get BET mask
        spm_orthviews('reposition', [0 0 0])
        aas_log(aap,false,'');
        
        % Mask current structural with the more accurate mask, to improve
        % future normalisation
        sV = spm_vol(Simg);
        sY = spm_read_vols(sV);
        sY = sY .* mY;
        spm_write_vol(sV,sY);
        
        %% DIAGNOSTIC
        subjname = aas_prepare_diagnostic(aap,subj);
        
        % This will only work for 1-7 segmentations
        OVERcolours = aas_colours;
        
        %% Draw native template
        spm_check_registration(Simg)
        
        % Add segmentations...
        for t = 1:(size(outSeg,1))
            spm_orthviews('addcolouredimage',1,outSeg(t,:), OVERcolours{t})
        end
        
        spm_orthviews('reposition', [0 0 0])
        
        print('-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' subjname '.jpeg']));
        
        %% Diagnostic VIDEO
        if aap.tasklist.currenttask.settings.diagnostic
            Ydims = {'X', 'Y', 'Z'};
            
            for d = 1:length(Ydims)
                aas_image_avi( Simg, ...
                    outSeg, ...
                    fullfile(aap.acq_details.root, 'diagnostics', [mfilename '__' subjname '_' Ydims{d} '.avi']), ...
                    d, ... % Axis
                    [800 600], ...
                    2, ... % Rotations
                    'none'); % No outline
            end
            try close(2); catch; end
        end
        
        % Another diagnostic image, looking at how well the segmentation  worked...
        Pthresh = 0.95;
        
        ROIdata = roi2hist(Simg, ...
            outSeg, Pthresh);
        
        [h, pv, ci, stats] = ttest2(ROIdata{1}, ROIdata{2});
        
        title(sprintf('GM vs WM... T-val: %0.2f (df = %d)', stats.tstat, stats.df))
        
        print('-djpeg','-r200',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' subjname '_Hist.jpeg']));
        
        % Now put our BETmask in the BETmask stream, but without deleting
        % the original BET mask...
        if aap.tasklist.currenttask.settings.brainMask
            aap=aas_desc_outputs(aap,subj,'BETmask', ...
                strvcat(V.fname, aas_getfiles_bystream(aap,subj,'BETmask')));
        else
            delete(BETmask)
        end
        aap=aas_desc_outputs(aap,subj,'structural',Simg);
        aap=aas_desc_outputs(aap,subj,'segmentation',outSeg);
end