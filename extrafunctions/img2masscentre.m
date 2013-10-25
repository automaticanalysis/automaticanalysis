function coords = img2masscentre(Mimg, binarise)
if nargin == 1
    binarise = 0;
end


% Load image first
V = spm_vol(Mimg);
M = spm_read_vols(V);

if binarise
    % Round the image...
    M = round(M);
end

sM = sum(M(:));

[X, Y, Z] = meshgrid(1:V.dim(1), 1:V.dim(2), 1:V.dim(3));

coords = ...
    [sum(X(:).*M(:))./sM, ...
    sum(Y(:).*M(:))./sM, ...
    sum(Z(:).*M(:))./sM];