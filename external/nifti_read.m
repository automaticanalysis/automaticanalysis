% Reads NIFTI files up to 5D

function Y = nifti_read(fname)
N = nifti(fname);
dim = [N.dat.dim 1]; inds = get_inds(dim(4:end));
for i = 1:size(inds,1)
    Y(:,:,:,i) = read3D(fname,inds(i,:));
end
Y = reshape(Y,dim);
end

function Y3D = read3D(fname,indices)
for i = numel(indices):-1:1
    fname = [fname ',' num2str(indices(i))];
end
Y3D = spm_read_vols(spm_vol(fname));
end

function inds = get_inds(dim)
[inds(:,3), inds(:,2), inds(:,1)] = ind2sub(dim,1:prod(dim));
for i = size(inds,2):-1:1
    if ~std(inds(:,i)), break; end
end
inds = inds(:,i+1:end);
if isempty(inds), inds = 1; end
end