function conv = FS_mri_convert(H,aap,varargin)
H = cell2mat(H);
[p, ia] = unique(spm_file({H.Filename},'path'));
fname = 'FSconv-.nii';
aas_runFScommand(aap,sprintf('mri_convert %s %s --split',H(ia(1)).Filename,fname));
conv.files = cellstr(spm_select('FPList',p{1},'^FSconv-.*'));