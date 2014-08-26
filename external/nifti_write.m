function nifti_write(fname,Y,desc,inf)
N      = nifti;
N.dat  = file_array(fname,inf.dim,inf.dt,0,inf.pinfo(1),inf.pinfo(2));
N.mat  = inf.mat;
N.mat0 = inf.mat;
N.mat_intent  = 'Scanner';
N.mat0_intent = 'Scanner';
N.descrip     = desc;
create(N);
N.dat(:,:,:) = Y;
end