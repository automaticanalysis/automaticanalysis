function cg_slice_overay_ui
%__________________________________________________________________________
% Christian Gaser
% $Id: cg_slice_overlay_ui.m 404 2011-04-11 10:03:40Z gaser $

OV.reference_image = fullfile(spm('dir'),'canonical','single_subj_T1.nii');
OV.reference_range = [0.05 0.6];                         % intensity range for reference image
OV.opacity = Inf;                                      % transparence value for overlay (<1)
OV.cmap    = jet;                                      % colormap for overlay

% name of files
OV.name=str2mat(fullfile(spm('dir'),'tpm','grey.nii'),...
                fullfile(spm('dir'),'tpm','white.nii'));
                
% range for each file
% Use range 0..0 if you want to autoscale range.
% If you are using log. scaling, check the highest p-value in the table
% and approximate the range; e.g. for a max. p-value of p=1e-7 and a
% threshold of p<0.001 use a range of [3 7]. Check cg_spmT2x.m for details.
% If you are unsure, simply use the autorange option by using a range of [0 0].
% The log-scaled values are calculated by -log10(1-p):
% p-value       -log10(1-P)
%  0.1           1
%  0.05          1.3
%  0.01          2
%  0.001         3
%  0.0001        4

% Number of fields in range should be the same as number of files (see above)
% or give one field, which is valid for all.
% Be carefule: intensities below the lower range are not shown!
OV.range   =[[0.5 1]; [0.5 1]];

% selection of slices and orientations
OV.slices_str = char('[-30 -5 5 40 50 60]','-30:2:30','-20:5:45');
OV.transform = char('axial','sagittal','coronal');
OV.printstr = 'print -r300 -painters -noui';

% define output format of slices
OV.labels.format = '%3.1f';

% Comment this out if you don't wish slice labels
%OV.labels = [];

% Comment this out if you don't wish colorbar
%OV.cbar = [];

cg_slice_overlay(OV)
