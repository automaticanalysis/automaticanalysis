function [aap,resp]=aamod_dartel_createtemplate(aap,task)
% AA module - Dartel create template
% Jonathan Peelle, Rik Henson, Rhodri Cusack
% MRC CBU Cambridge Oct 2008
%
% This script creates a template using DARTEL. The default behavior
% in DARTEL is to leave created templates in the directory of the
% first subject in the list of imported images.  This script copies
% these templates to templates_[studyname], where studyname is
% aap.directory_conventions.dartelsubjdirname.
%
% Template creation can be quite time consuming; for example, 48
% hours for a template of ~80 subjects is not unreasonable.
%
% See aamod_dartel_import or aamod_segment8 for more information.
%
% $Id$

resp='';

% possible tasks 'doit','report','checkrequirements'
switch task
 case 'domain'
  resp='study';
 case 'description'
  resp='Create DARTEL template';
 case 'doit'
  % report which version of spm_dartel_template we are using
  aas_log(aap, false, sprintf('Using %s.\n', which('spm_dartel_template')));
     
  % number of tissue classes included (usually 2 or 6)
  numtissues = aap.tasksettings.aamod_dartel_createtemplate.numtissueclasses; 
                                                                 
                                                        
  % keep track of any we couldn't find
  notfound = {};
  
  % get images
  clear allimages
  for i=1:length(aap.acq_details.subjects)
    dir_subj= aas_getsubjpath(aap,i);
    dir_subjdartel=fullfile(dir_subj,aap.directory_conventions.structdirname,aap.directory_conventions.dartelsubjdirname);
    for k=1:numtissues
      imgs=spm_select('fplist', dir_subjdartel, sprintf('^rc%d.*nii',k));
      
      % Note any subjects for which we couldn't find an image
      if isempty(imgs) || strcmp(imgs, '/')
        notfound = {notfound{:} aap.acq_details.subjects(i).mriname};
      end
      
      if (size(imgs,1)~=1)
        aas_log(aap,true,sprintf('Did not find exactly 1 rc%d image in %s',k,dir_subjdartel));
      end
      allimages{k}{i}=imgs;
    end
  end
  
  % Send an error if any required images not found
  if ~isempty(notfound)
    for m=1:length(notfound)
      fprintf('Did not find an image for %s.\n', notfound{m});
    end
    error('Did not find all required images.')
  end
  

  % Set up job
  % below based on tbx_cfg_dartel 15 July 2009
  % eventually these should probably be user-definable
  param = struct(...
      'its',{3,3,3,3,3,3},...
      'rparam',{[4 2 1e-6],[2 1 1e-6],[1 0.5 1e-6],...
                [0.5 0.25 1e-6],[0.25 0.125 1e-6],[0.25 0.125 1e-6]},...
      'K',{0,0,1,2,4,6},...
      'slam',{16,8,4,2,1,0.5});
  
  settings = struct('template','Template','rform',aap.tasklist.currenttask.settings.rform,...
                    'param',param,...
                    'optim', struct('lmreg',0.01, 'cyc', 3, 'its', 3));
  
  

  job = struct('images',{allimages}, 'settings',settings);

  % run the script
  spm_dartel_template(job);
        
  
  % Make a folder in the main study directory to hold templates
  templatedir = fullfile(aap.acq_details.root, sprintf('templates_%s',aap.directory_conventions.dartelsubjdirname));
  if ~isdir(templatedir)
    mkdir(templatedir);
  end
  
  
  dir_subj= aas_getsubjpath(aap,1);               
  dir_struct=fullfile(dir_subj,aap.directory_conventions.structdirname);
  dir_subjdartel=fullfile(dir_struct,aap.directory_conventions.dartelsubjdirname);
  
  try
    fprintf('Moving templates from %s to %s...', dir_subjdartel, templatedir);
    system(sprintf('mv %s/Template*.nii %s/', dir_subjdartel, templatedir));
    fprintf('done.\n')
  catch
    aas_log(aap,true,sprintf('Error moving templates to %s, but everything else should be ok.\n', templatedir))
  end
end









