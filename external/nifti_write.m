function nifti_write(fname,Y,desc,inf)
N      = nifti;
N.dat  = file_array(fname,inf.dim,inf.dt,0,0,inf.pinfo(2));
N.mat  = inf.mat;
N.mat0 = inf.mat;
N.mat_intent  = inf.private.mat_intent;
N.mat0_intent  = inf.private.mat0_intent;
N.descrip     = desc;
create(N);
N.dat(:,:,:) = Y;
end