% Gets voxel and matrix size in mm
function [mmVox mmFOV] = vox2mm(V)

mmVox = zeros(1,3);

% Get the matrix
M = V.mat;

% Get the field of view in voxels
FOV = V.dim;

for d = 1:3
    dM = [0 0 0; 0 0 0];
    dM(1,d) = 1;
    
    pM = cor2mni(dM, M);
    mmVox(d) = pdist(pM);    
end

mmFOV = FOV .* mmVox;

end