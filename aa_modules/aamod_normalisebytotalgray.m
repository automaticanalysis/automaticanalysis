% AA module - Divide voxelwise GM by TGM as proportional scaling for VBM
% [aap,resp]=aamod_normalisebytotalgray(aap,task,i)
% Rhodri Cusack MRC CBU Cambridge Oct 2008
% Divide structural images by total gray matter volume as a form of
% proportinal scaling. Output images have a 'g' prepended.
% Will work on segmented structural images, or GM images processed
% through the DARTEL AA pipeline.
% Requires aamod_structuralstats to be run first to get total gray
% matter.
% $Id: aamod_normalisebytotalgray.m 159 2009-08-18 13:15:14Z jpeelle $
% @@@ THIS IS NOT YET TRANSFORMED TO AA4 @@@

function [aap,resp]=aamod_normalisebytotalgray(aap,task,i)

resp='';

% possible tasks 'doit','report','checkrequirements'
switch task
 case 'domain'
  resp='subject';
 case 'report'
  resp='Noramlise images by total gray matter.';
  
 case 'doit'
  dir_subj= aas_getsubjpath(aap,i);
  dir_struct=fullfile(dir_subj,aap.directory_conventions.structdirname);
  load(fullfile(dir_struct, 'structuralstats.mat'));
  tgv=S.parts.mm3(1);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Divide segmented structural image
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  filters={fullfile(dir_struct,'mwc1*nii'), fullfile(dir_struct,'smwc1*nii')};

  for f=1:length(filters)
      fn=dir(filters{f});
      if length(fn)>0
          V=spm_vol(fullfile(dir_struct,fn.name));
          Y=spm_read_vols(V);
          Y=Y/tgv;
          [pth fle ext]=fileparts(V.fname);
          V.fname=fullfile(pth,['g' fle ext]);
          spm_write_vol(V,Y);
      end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % If DARTEL directory exists, also divide that
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if isfield(aap.directory_conventions, 'dartelsubjdirname') && ~isempty(aap.directory_conventions.dartelsubjdirname)
    dir_dartel = fullfile(dir_struct, aap.directory_conventions.dartelsubjdirname);
    filters = {fullfile(dir_dartel,'mwrc1*nii'), fullfile(dir_dartel, 'smwrc1*nii')};
    
    for f=1:length(filters)
        fn = dir(filters{f});
        if length(fn)>0
            V=spm_vol(fullfile(dir_dartel,fn.name));
            Y=spm_read_vols(V);
            Y=Y/tgv;
            [pth fle ext]=fileparts(V.fname);
            V.fname=fullfile(pth,['g' fle ext]);
            spm_write_vol(V,Y);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % If mni* files exist, also divide those
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mnifiles = spm_select('fplist', dir_dartel, '^mnimwrc1');
    if ~isempty(mnifiles) && ~strcmp(mnifiles, '/')
      % we found a file
      if size(mnifiles,1)~=1
        aas_log(aap, true, sprintf('Error looking for MNI-normalised DARTEL GM file'))
      end
      
      V = spm_vol(mnifiles);
      Y = spm_read_vols(V);
      Y = Y/tgv;
      
      [pth fle ext]=fileparts(V.fname);
      V.fname=fullfile(pth,['g' fle ext]);
      spm_write_vol(V,Y);
      
    end % checking for MNI files
  end % DARTEL scaling
end