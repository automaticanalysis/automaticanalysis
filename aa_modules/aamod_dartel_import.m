function [aap,resp]=aamod_dartel_import(aap,task,i)
% AA module - Dartel import
% Jonathan Peelle, Rik Henson, Rhodri Cusack
% MRC CBU Cambridge Oct 2008
%
% The DARTEL AA modules were designed to essentially follow John
% Ashburner's DARTEL guide:
%
% http://www.fil.ion.ucl.ac.uk/~john/misc/dartel_guide.pdf
%
% In general these steps are:
%  1) segment images into GM and WM (aamod_norm_noss)
%  2) import these segmentations into DARGEL (aamod_dartel_import)
%  3) Create a template (aamod_dartel_createtemplate)
%  4) Write normalized images (aamod_dartel_normwrite)
%
%  5) (optional) calculate template-to-MNI transform (aamod_dartel_mnitransform)
%
%  6) (optional) write images to MNI space
%  (aamod_dartel_all2mni)
%
%  7) (optional) Proportionally scale images by total GM
%  (aamod_normalisebytotalgray) (requires aamod_structuralstats to
%  have been run previously)
%
%  8) Smooth structurals (aamod_smoothstructurals)
%
%
% Because DARTEL creates a custom template it has to be re-run for
% each group of subjects. Each time you create a template, a
% subdirectory is created in each subject's structurals directory,
% based on aap.directory_conventions.dartelsubjdirname.
%
%
% $Id$

resp='';

% possible tasks 'doit','report','checkrequirements'
switch task
 case 'domain'
  resp='subject';
 case 'description'
  resp='Import images for DARTEL';
  
 case 'doit'
  dir_subj= aas_getsubjpath(aap,i);
  dir_struct=fullfile(dir_subj,aap.directory_conventions.structdirname);
  dir_subjdartel=fullfile(dir_struct,aap.directory_conventions.dartelsubjdirname);
  aas_makedir(aap,dir_subjdartel);
  filter_snmat=fullfile(dir_struct,'*seg_sn.mat');
  
  fn_snmat=dir(filter_snmat);
  if (length(fn_snmat)<1)
    aas_log(aap,true,sprintf('failed to find snmat with filter %s',snmat_filter));
  end
  
  clear initial;
  initial.matnames = {fullfile(dir_struct,fn_snmat(1).name)};
  initial.odir = cellstr(dir_subjdartel);
  initial.bb = nan(2,3);
  initial.vox = aap.tasklist.currenttask.settings.vox;
  initial.image = 0;
  initial.GM = 1;
  initial.WM = 1;
  initial.CSF = 0;
  
  spm_dartel_import(initial);
  
end









