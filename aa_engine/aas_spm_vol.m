% Like spm_vol, except that it will read .nii.gz files
% Useful for FSL integration
% Rhodri Cusack, London, Canada, June 2013

function [V Y XYZ]=aas_spm_vol(fn,dontdelete)
if ~exist('dontdelete','var')
    dontdelete=false;
end;

[pth nme ext]=aas_fileparts(fn);

if strcmp(ext,'.nii.gz')
    % Uncompress and read
    [temp_pth temp_nme temp_ext]=fileparts(tempname);
    tf=fullfile(pth,[temp_nme ext]);
    copyfile(deblank(fn),tf);
    [s w]=aas_shell(['gunzip ' tf]);
    V=spm_vol(tf(1:end-3));
    if nargout==3
        [Y XYZ]=spm_read_vols(V);
    else
        Y=spm_read_vols(V);
    end;
    if ~dontdelete
        delete(tf(1:end-3));
    end;
else
    V=spm_vol(fn);
    if nargout==3
        [Y XYZ]=spm_read_vols(V);
    else
        Y=spm_read_vols(V);
    end;
end;