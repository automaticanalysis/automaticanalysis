function [ext]= aas_getfslext(aap)

switch (aap.directory_conventions.fsloutputtype)
    case 'NIFTI'
        ext='.nii';
    case 'NIFTI_GZ'
        ext='.nii.gz';
end;
end