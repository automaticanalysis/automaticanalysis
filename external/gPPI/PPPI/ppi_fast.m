function ppi_fast(fsfastparams,ppiparams)
%% This function will port fsfast structure to SPM.mat and then run the PPI process in SPM using PPPI
% 
% For details on ppiparams, please see PPPI.m
%
%  
% Written by Donald McLaren (mclaren@nmr.mgh.harvard.edu)
% GRECC, Edith Norse Roers Memorial Veterans Hospital, Bedford, MA
% Department of Neurology, Massachusetts General Hospital and Harvard
%   Medical School
% 1/30/2011

%% Check SPM version
try
    if ~strcmp(spm('Ver'),'SPM8')
    disp('PROGRAM ABORTED:')
    disp('  You must use SPM8 to process your data; however, you can use SPM.mat files')
    disp('  generated with SPM2 or SPM5. In these cases, simply specify the option SPMver')
    disp('  in single qoutes followed by a comma and the version number.')
    disp(' ')
    disp('Make sure to add SPM8 to your MATLAB path before re-running.')
    return
    end
catch
    disp('SPM not initialized');
    return
end

try
    [SPM,fsfastparams]=spm_fsfastfirstlevel_checkparams(fsfastparams);
catch
    error('Parameter Structure has wrong values.');
end

%Check for Contrasts
% if ~isfield(ppiparams,'Contrasts') || isempty(ppiparams.Contrasts)
%        for ii=1:numel(fsfastparams.runflac(1).flac.con)
%            if any(strmp(fsfastparams.runflac(1).flac.con(ii).name,fsfastparams.Contrasts(ii).name))
%                ind=find(strmp(fsfastparams.runflac(1).flac.con(ii).name,fsfastparams.Contrasts(ii).name));
%            else
%                ind=numel(fsfastparams.Contrasts)+1;
%            end
%            ppiparams.Contrasts(ind).name=fsfastparams.runflac(1).flac.con(ii).name;
%            ppiparams.Contrasts(ind).STAT='T';
%            ppiparams.Contrasts(ind).Weighted=0;
%            if ~isfield(ppiparams.Contrasts(ind),'left'); ppiparams.Contrasts(ind).left={}; end
%            if ~isfield(ppiparams.Contrasts(ind),'right'); ppiparams.Contrasts(ind).right={}; end
%            tmp=find(fsfastparams.runflac(jj).flac.con(ii).c~=0);
%            for kk=1:numel(tmp)
%                if tmp(kk)>0
%                    ppiparams.Contrasts(ind).left{end+1}=SPM.Sess(jj).U(tmp(kk)).name{1};
%                else
%                    ppiparams.Contrasts(ind).right{end+1}=SPM.Sess(jj).U(tmp(kk)).name{1};
%                end
%            end
%        end
%        for ii=1:numel(ppiparams.Contrasts)
%            ppiparams.Contrasts(ii).left=unique(ppiparams.Contrasts(ii).left);
%            ppiparams.Contrasts(ii).right=unique(ppiparams.Contrasts(ii).right);
%            if isempty(ppiparams.Contrasts(ii).left); ppiparams.Contrasts(ii).left={'none'}; end
%            if isempty(ppiparams.Contrasts(ii).right); ppiparams.Contrasts(ii).right={'none'}; end
%        end
% end

%% Set Directories
tt = pwd;
try
    cd(SPM.swd);
catch
    mkdir(SPM.swd);
    cd(SPM.swd);
end 
if fsfastparams.estimate==1
     disp('Estimating Task First Level Model.') 
     try delete beta_00*; end
     try delete ResMS.*; end
     try delete RPV.*; end
     try delete mask.*; end
     
     %Build Design Matrix
     SPM = spm_fmri_spm_ui(SPM); % Build Design Matrix
     SPM.xM.T(:) = -Inf;  %% disable threshold masking
     try
         SPM.xM.VM = spm_vol(fsfastparams.runflac(1).flac.maskfspec);
     catch
         %error('No mask file.');
     end
     save SPM SPM
     SPM = spm_spm(SPM); % Estimate Model
     %try
     %   SPM=spm_firstlevel_contrasts(ppiparams)
     %catch
     %end
     disp('Task First Level Model has now been estimated. Progam will now begin PPI.')
else
     disp('Task First Level Model already estimated. Progam will now begin PPI.')
end
PPPI(ppiparams);

