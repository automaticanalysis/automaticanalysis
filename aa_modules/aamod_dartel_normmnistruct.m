function [aap,resp]=aamod_dartel_normmnistruct(aap,task)
% After DARTEL template has been created, write out normalized
% versions (MNI space) of each subject's structural image.
%
% This function will also do smoothing, as specified in
% aap.tasksettings.aamod_dartel_normmnistruct.fwhm. The default
% (specified in the .xml file) is .5 mm (to help with aliasing).
%
% See other aamod_dartel functions for more information.
%
% Jonathan Peelle
% $Id$


resp='';

% possible tasks 'doit','report','checkrequirements'
switch task
 case 'domain'
  resp='study';
 case 'report'
  resp='Write smoothed normalised structural images in MNI space for DARTEL.'
 case 'doit'

  % template
  templatedir = fullfile(aap.acq_details.root, ['templates_' aap.directory_conventions.dartelsubjdirname]);
  template = spm_select('fplist', templatedir, '^Template_6\.nii$');
  if isempty(template) || strcmp(template, '/')
    aas_log(aap,true,'Did not find Template6.nii.');
  end
  
  % initialize images
  allimages = {};
  
  for i=1:length(aap.acq_details.subjects)
    dir_subj= aas_getsubjpath(aap,i);
    dir_struct=fullfile(dir_subj,aap.directory_conventions.structdirname);
    dir_subjdartel=fullfile(dir_struct,aap.directory_conventions.dartelsubjdirname);
    
    % flow fields..
    job.data.subj(i).flowfield{1} = spm_select('fplist', dir_subjdartel, '^u_.*nii');
    
    % get the structural image to normalize, if requested
    imgs = spm_select('fplist', dir_struct, sprintf('^%s.*nii', aap.tasksettings.aamod_dartel_normmnistruct.structprefix));

      job.data.subj(i).images = cellstr(imgs);
  end % going through subjects
  
  % set up job
  job.template{1} = template;
  job.bb = nan(2,3);
  job.vox = ones(1,3) * aap.tasksettings.aamod_dartel_normmnistruct.vox;
  job.fwhm = aap.tasksettings.aamod_dartel_normmnistruct.fwhm;
  job.preserve = aap.tasksettings.aamod_dartel_normmnistruct.preserve;
  
  aas_log(aap, false, sprintf('Running with %s...', which('spm_dartel_norm_fun')));
  spm_dartel_norm_fun(job);
    
  
  % should we add diagnostic image writing at some point?
end
