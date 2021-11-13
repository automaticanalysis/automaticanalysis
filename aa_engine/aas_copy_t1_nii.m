function tmpfile = aas_copy_t1_nii(aap, localpath)
% Copy a (unzipped) T1 template to localpath
%
% First try FSL's MNI152_T1_1mm, fallback is SPM provided T1template
%
% INPUTS
% - aap: aa parameters structure
% - localpath: directory to put the T1 template file
%
% OUTPUTS
% - tmpfile: full path to the local copy of the T1 template

tmpfile = fullfile(aap.directory_conventions.fsldir,'data','standard','MNI152_T1_1mm.nii.gz'); % use FSL highres
if exist(tmpfile,'file')
    gunzip(tmpfile, localpath);
    tmpfile = fullfile(localpath, 'MNI152_T1_1mm.nii');
else % use SPM's
    [~, SPMtool] = aas_cache_get(aap,'spm');
    tmpfile = fullfile(SPMtool.toolPath, aap.directory_conventions.T1template);
    copyfile(tmpfile, localpath);
    tmpfile = spm_file(tmpfile, 'path', localpath);
end

end