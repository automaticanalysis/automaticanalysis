function out = cg_vbm8_longitudinal_multi_run(job)
% Call cg_vbm8_longitudinal for multiple subjects
%
% Christian Gaser
% $Id: cg_vbm8_longitudinal_multi_run.m 431 2011-09-30 14:18:03Z gaser $

global data_long

warning off;

for i=1:numel(job.subj),
    out(i).files = cell(numel(job.subj(i).mov),1);
    m = numel(job.subj(i).mov);
    data = cell(m,1);
    for j=1:m
        [pth,nam,ext,num] = spm_fileparts(job.subj(i).mov{j});
        out(i).files{j} = fullfile(pth,['wp1mr', nam, ext, num]);
        data{j} = job.subj(i).mov{j};
    end
    data_long = {data};
    cg_vbm8_longitudinal;
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
end;

warning on;
