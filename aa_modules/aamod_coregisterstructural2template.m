function [aap,resp]=aamod_coregisterstructural2template(aap,task,subjind)
% AAMOD_COREGISTERSTRUCTURAL2TEMPLATE Coregister structural image to
% template.
%
% This can aid segmentation, for example, by improving starting estimates.
%

resp='';

% possible tasks 'doit','report','checkrequirements'
switch task
  case 'domain'
    resp='subject';
  case 'report'
    resp='Coregister structural image with T1 template.'
  case 'doit'
    
    % get the template
    template = fullfile(spm('dir'), 'canonical', 'avg152T1.nii');
    if ~exist(template)
      aas_log(aap, true, sprintf('Couldn''t find template T1 image %s.', template));
    end
    
    
    % get the structural image
    structImg = aas_getfiles_bystream(aap, subjind, 'structural');
        
    
    VG = spm_vol(template);
    VF = spm_vol(structImg);
    
    x = spm_coreg(VG, VF);
    M = inv(spm_matrix(x));
    
    PO = structImg;
    
    % TODO: add option to move functional images along with
    %     if cfg.move_functional
    %       PO = strvcat(PO, funimages);
    %     end
    
    MM = zeros(4,4,size(PO,1));
    
    for j=1:size(PO,1)
      MM(:,:,j) = spm_get_space(deblank(PO(j,:)));
    end
    
    for j=1:size(PO,1)
      spm_get_space(deblank(PO(j,:)), M*MM(:,:,j));
    end
   
    % describe outputs
    aap = aas_desc_outputs(aap, subjind, 'structural', structImg);
    
end

