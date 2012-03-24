% AA module - smoothing for structural images
% function [aap,resp]=aamod_smooth_structurals(aap,task,i)
% Kernel size determined by aap.spm_analysis.FWHM (default 12mm).
% Rhodri Cusack MRC CBU Cambridge October 2008
% $Id: aamod_smooth_structurals.m 159 2009-08-18 13:15:14Z jpeelle $

function [aap,resp]=aamod_smooth_structurals(aap,task,i)

resp='';

switch task
    case 'domain'
        resp='subject';

    case 'description'
        resp='Smooth structural images.';

    case 'summary'
        subjpath=aas_getsubjpath(i);
        resp=sprintf('Smooth run on %s\n',subjpath);
        
    case 'report'
     
    case 'doit'

        subj_dir = aas_getsubjpath(aap,i);
        
        darteldir = aap.directory_conventions.dartelsubjdirname;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Images to smooth
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % (basically try to smooth any structurals found:
        % segmentations, normalized segmentations, normalized
        % dartel images, globally-scaled images, or MNI-scaled
        % DARTEL images (mni*).

        % saved space by not specifying tissue type - wc* instead of wc1*
        % and wc2*.
        filters={'wc*.nii','mwc*.nii','g*.nii', fullfile(darteldir,'mwrc*.nii'),fullfile(darteldir,'wrc*nii'), fullfile(darteldir,'mni*nii'), fullfile(darteldir,'g*.nii')};        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get images
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        imgs = [];
        for f=1:length(filters)
          P = dir(fullfile(subj_dir,aap.directory_conventions.structdirname,filters{f}));
          [filterpth fle ext]=fileparts(filters{f});
          for b=1:length(P)
            imgs = strvcat(imgs, fullfile(subj_dir,aap.directory_conventions.structdirname,filterpth,P(b).name));
          end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Smooth
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        s   = aap.tasklist.currenttask.settings.FWHM;
        for f = 1:size(imgs,1)
          [pth,nam,xt] = spm_fileparts(imgs(f,:));
          U = fullfile(pth,['s' nam xt]);
          spm_smooth(imgs(f,:),U,s);
        end
        
    case 'checkrequirements'

    otherwise
     aas_log(aap,1,sprintf('Unknown task %s',task));
end