% AA module - Structural statistics for VBM
% Rhodri Cusack MRC CBU Cambridge Oct 2008
% Get volume in mm^3 of GM, WM, CSF, and total intracranial volume (TIV).
% The default prefix is mwc* (e.g., volume for mwc1*, mwc2*, and mwc3*
% images).  However, this can be changed in:
%  aap.tasksettings.aamod_structuralstats.prefix
% if you want this to be something else.
% $Id: aamod_structuralstats.m 165 2009-11-12 16:54:14Z jpeelle $

function [aap,resp]=aamod_structuralstats(aap,task,i)
resp='';

% possible tasks 'doit','report','checkrequirements'
switch task
 case 'domain'
  resp='subject';
 case 'report'
  resp='Get statistics on segmented structural images.'  
 case 'doit'
  dir_subj = aas_getsubjpath(aap,i);
  dir_struct = fullfile(dir_subj,aap.directory_conventions.structdirname);
  
  %dir_seg8 = fullfile(dir_struct, 'seg8');
  %isseg8 = isdir(dir_seg8); % only try this if seg8 directory exists
  
  parts={'1=grey matter','2=white matter','3=csf'};

  for k=1:3
    filter_part=fullfile(dir_struct,sprintf('%s%d*',aap.tasksettings.aamod_structuralstats.prefix,k));
    dir_part=dir(filter_part);
    
    % If we can't find a file, that's ok as long as something is in
    % the segment8 directory
    if (length(dir_part)~=1)
      aas_log(aap,true,sprintf('Error looking for modulated normalised file %s',filter_part));
    else
      V=spm_vol(fullfile(dir_struct,dir_part.name));
      Y=spm_read_vols(V);
      spacedesc=spm_imatrix(V.mat);
      volume=prod(abs(spacedesc(7:9)));
      S.parts.desc=parts;
      Y=Y(~isnan(Y));
      S.parts.vox(k)=sum(Y(:));
      S.parts.mm3(k)=S.parts.vox(k)*volume;
    end
    
%     % stats for segmentation8
%     if isseg8
%       filter_part_seg8=fullfile(dir_seg8,sprintf('%s%d*',aap.tasksettings.aamod_structuralstats.prefix,k));
%       dir_part_seg8=dir(filter_part_seg8);
%       
%       V=spm_vol(fullfile(dir_seg8,dir_part_seg8.name));
%       Y=spm_read_vols(V);
%       spacedesc=spm_imatrix(V.mat);
%       volume=prod(abs(spacedesc(7:9)));
%       S_seg8.parts.desc=parts;
%       Y=Y(~isnan(Y));
%       S_seg8.parts.vox(k)=sum(Y(:));
%       S_seg8.parts.mm3(k)=S_seg8.parts.vox(k)*volume;
%     end 
    
  end % going through tissue classes

  S.TIV.mm3=sum(S.parts.mm3);
  S.TIV.vox=sum(S.parts.vox);
  save(fullfile(dir_struct,'structuralstats.mat'),'S');
  
%   if isseg8
%     S_seg8.TIV.mm3=sum(S_seg8.parts.mm3);
%     S_seg8.TIV.vox=sum(S_seg8.parts.vox);
%     save(fullfile(dir_seg8,'structuralstats_seg8.mat'),'S_seg8');    
%   end
end









