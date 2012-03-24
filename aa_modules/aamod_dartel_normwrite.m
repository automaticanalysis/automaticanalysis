function [aap,resp]=aamod_dartel_normwrite(aap,task,i)
% AA module - Dartel import
% Jonathan Peelle, Rik Henson, Rhodri Cusack
% MRC CBU Cambridge Oct 2008
%
% After DARTEL template has been created, write out normalized
% versions of each subject's images.
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
  resp='Write normalised images for DARTEL.'
 case 'doit'
  % directories...
  dir_subj= aas_getsubjpath(aap,i);
  dir_struct=fullfile(dir_subj,aap.directory_conventions.structdirname);
  dir_subjdartel=fullfile(dir_struct,aap.directory_conventions.dartelsubjdirname);
  
  % flow fields..
  flowfieldimages = cellstr(spm_select('fplist', dir_subjdartel, '^u_'));
  for k=1:length(flowfieldimages)
    flowfieldimages{k} = strtok(flowfieldimages{k},',');
  end
  
  % images...
  clear allimages
  for k=1:2
    imgs=spm_select('fplist', dir_subjdartel, sprintf('^rc%d.*nii',k));
    if (size(imgs,1)~=1)
      aas_log(aap,true,sprintf('Did not find exactly 1 rc%d image in %s',k,dir_subjdartel));
    end
    allimages{1}{k}=imgs;
  end
  
  % set up job
  job.flowfields = flowfieldimages;
  job.images = allimages;
  job.jactransf = 1;
  job.K = 6;
  job.interp = 1;
  
  % spm_dartel_norm saves images to the current working directory
  originaldir = pwd;
  cd(dir_subjdartel);
  spm_dartel_norm(job);
  cd(originaldir);
  
  % should add diagnostic image writing?
end
