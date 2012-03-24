function [aap,resp]=aamod_dartel_all2mni(aap,task,i)
% AA module - Dartel conversion to MNI space (modulated)
% Jonathan Peelle
% MRC CBU Cambridge January 2009
%
% Based on John Ashburner's email:
%
% https://www.jiscmail.ac.uk/cgi-bin/wa.exe?A2=ind0803&L=SPM&P=R22666 
%
% See aamod_dartel_import for more information.
%
% $Id$

resp='';

% possible tasks 'doit','report','checkrequirements'
switch task
 case 'domain'
  resp='subject';
 case 'report'
  resp='Warp DARTEL images to MNI space, modulating to preserve total GM.'
 case 'doit'
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Get the normalization parameters
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  param_file = fullfile(aap.acq_details.root, sprintf('templates_%s', aap.directory_conventions.dartelsubjdirname), 'Template_6_sn.mat');
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Get the image for this subject
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  dir_subj= aas_getsubjpath(aap,i);
  dir_struct=fullfile(dir_subj,aap.directory_conventions.structdirname);
  dir_subjdartel=fullfile(dir_struct,aap.directory_conventions.dartelsubjdirname);
  
  % select the original registered segmentation
  dartelname = spm_select('fplist', dir_subjdartel, '^mwrc1');
  
  if isempty(dartelname) || strcmp(dartelname, '/')
    aas_log(aap, true, 'Problem finding the dartel image.');
  end
  
  
  % Copy this image so we have the original, as well as the new
  % header-tweaked version
  [pth, nm, ext] = fileparts(dartelname);
  
  newdartel = fullfile(pth, ['mni' nm ext]);
  system(sprintf('cp %s %s', dartelname, newdartel));
  dartelname = newdartel;
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Normalize and modulate by tweaking header info
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % (This is from code below, sent by John Ashburner - JP)
  
  % PN = spm_select(1,'.*_sn.mat','Select sn.mat file');
  % (we already have this in param_file)
  
  % PI = spm_select(inf,'nifti','Select images');
  % (we already have this as dartelname)
  
  % % Determine affine transform from header
  sn    = load(deblank(param_file));
  M     = sn.VG(1).mat/(sn.VF(1).mat*sn.Affine);
  
  
  % Scaling by inverse of Jacobian determinant, so that
  % total tissue volumes are preserved.
  scale = 1/abs(det(M(1:3,1:3)));
  
  
  % Pre-multiply existing headers by affine transform
  Ni = nifti(deblank(dartelname));
  
  % Pre-multiply existing header by affine transform
  Ni.mat = M*Ni.mat;
  Ni.mat_intent='MNI152';
  
  % Change the scalefactor.  This is like doing a "modulation"
  Ni.dat.scl_slope = Ni.dat.scl_slope*scale;
  
  % Write the file
  create(Ni);                          
end

