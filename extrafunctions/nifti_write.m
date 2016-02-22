% Write NIFTI files up to 5D

function nifti_write(fname,Y,desc,V)
dim = size(Y);
N      = nifti;
N.dat  = file_array(fname,dim,V.dt);
N.mat  = V.private.mat;
N.mat0 = V.private.mat;
N.descrip     = desc;
create(N);

dim = [dim 1];
for i = 1:prod(dim(4:end))
    N.dat(:,:,:,i) = Y(:,:,:,i);   
    spm_get_space([N.dat.fname ',' num2str(i)], V.mat);
end
N.dat = reshape(N.dat,dim);
end