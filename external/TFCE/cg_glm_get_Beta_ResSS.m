function [Beta, ResMS] = cg_glm_get_Beta_ResSS(V,Vmask,X,pKX,TH,W);

fprintf('Compiling cg_glm_get_Beta_ResSS.c\n');

spm_dir = spm('dir');
if strcmp(spm('ver'),'SPM2')
	spm_src = spm_dir;
else
	spm_src = fullfile(spm_dir,'src');
end

c_dir = which('cg_glm_get_Beta_ResSS.c');
disp(['mex -I' spm_src ' -outdir ' fileparts(c_dir) ' ' c_dir ' ' spm_src filesep 'spm_vol_utils.' mexext '.a'])
eval(['mex -I' spm_src ' -outdir ' fileparts(c_dir) ' ' c_dir ' ' spm_src filesep 'spm_vol_utils.' mexext '.a'])

try
	[Beta, ResMS] = cg_glm_get_Beta_ResSS(V,Vmask,X,pKX,TH,W);
end