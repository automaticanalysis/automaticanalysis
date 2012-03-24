function [aap,resp]=aamod_dartel_mnitransform(aap,task)
% AA module for determining MNI transformation for DARTEL template
%
% Jonathan Peelle
% MRC CBU January 2009
% based on 1.2.1 of John Ashburner's DARTEL guide
%
% AA DARTEL scripts assume template has been put in its own
% directory using aamod_dartel_createtemplate.  (The DARTEL default
% is to leave the template in the directory of the first subject.)
%
% See aamod_dartel_import for more information.
%
% $Id$

resp='';

% possible tasks 'doit','report','checkrequirements'
switch task
 case 'domain'
  resp='study';
 case 'report'
  resp='Calcualte MNI transform for DARTEL template.'
  
 case 'doit'
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Get the DARTEL template and GM template supplied by SPM
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %         % The template should be in the first subject's directory
  %         dir_subj= aas_getsubjpath(aap,1);
  %
  %         dir_struct=fullfile(dir_subj,aap.directory_conventions.structdirname);
  %         templatedir=fullfile(dir_struct,aap.directory_conventions.dartelsubjdirname);
  
  
  % new versions of aa for dartel put templates in their own
  % directory
  templatedir = fullfile(aap.acq_details.root, sprintf('templates_%s',aap.directory_conventions.dartelsubjdirname));
  
  % Get the final template from DARTEL - assumes it will be
  % Template_6.                        
  dartel_template = fullfile(templatedir, 'Template_6.nii');
  if ~exist(dartel_template)
    aas_log(aap,true,sprintf('Failed to find DARTEL template here: %s', dartel_template));
  end
  
  % Get the GM TPM from spm
  gm_template = fullfile(spm('Dir'), 'tpm', 'grey.nii');
  if ~exist(gm_template)
    aas_log(aap,true,sprintf('Failed to find GM TPM here: %s', gm_template));
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
  % Normalize this to the GM priors (Normalise: Estimate)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Note: SPM_VOL returned a 2 element array for the dartel template,
  % which messed things up.  SPM_VOL_NIFTI seems to work fine, so I'm
  % not sure what the issue was. Generally the DARTEL template
  % should have 2 volumes, one for GM and one for WM.
  
  VG = spm_vol_nifti(gm_template);
  VF = spm_vol_nifti(dartel_template);
  
  
  % The output file for the normalization routine
  [pth, nm] = fileparts(dartel_template);
  norm_file = fullfile(pth, [nm '_sn.mat']);
  
  flags = struct();
  flags.smosrc = 8;       % smooth source at 8 mm to match the GM template TPM
  flags.smoref = 0;       % don't smooth tpm as this is already smoothed
  flags.nits = 0;         % 0 nonlinear iterations = affine-only spatial normalisation
  flags.regtype = 'mni'; 
  
  params = spm_normalise(VG, VF, norm_file, [], [], flags);

end









