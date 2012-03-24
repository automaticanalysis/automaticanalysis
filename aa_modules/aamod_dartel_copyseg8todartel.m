function [aap,resp] = aamod_dartel_copyseg8todartel(aap,task,i)
%AAMOD_COPYSEGM8TODARTEL Bridges SEGMENT8 to DARTEL.
%
% In 'standard' DARTEL analysis, an 'import' step takes segmented images
% (c*) and turns them into rc* images; in AA this is done separately for
% each DARTEL analysis, and thus the rc* images are put in a separate
% directory each time.
%
% When segment8 is used (aamod_segment8) for updated segmentation routines,
% rc* images are automatically generated, saved in the seg8 folder in each
% subject's structurals folder. This module creates softlinks from each new
% DARTEL analysis folder to the rc* images so that you don't have to
% run segmentation8 for each new DARTEL analysis.

% Jonathan Peelle
% MRC Cognition and Brain Sciences Unit

resp='';

switch task
 case 'domain'
  resp='subject';
 case 'description'
  resp='SPM8 softlinks to seg8/rc* images from DARTEL directory.'
 case 'doit'
  
  global defaults;

  
  % structural directory for this subject
  subjdir = aas_getsubjpath(aap,i);
  structdir = fullfile(subjdir, aap.directory_conventions.structdirname);

  % DARTEL directory
  darteldir = fullfile(structdir, aap.directory_conventions.dartelsubjdirname);
  
  if ~isdir(darteldir)
      mkdir(darteldir);
  end
  
  
  % get all the rc* images from the seg8dir
  rcimages = spm_select('fplist', structdir, '^rc.*\.nii');
  
  % error if none found
  if size(rcimages,1) == 0
    aas_log(aap, 1, 'No rc* images found.');
  end
  
  for i=1:size(rcimages,1)
     img = strtok(rcimages(i,:)); % get this image       
     fprintf('Creating softlink to %s...', img);
     [pth, nm, ext] = fileparts(img);
     system(sprintf('ln -s %s %s', img, fullfile(darteldir,[nm ext])));
     fprintf('done.\n');
  end
  
end

