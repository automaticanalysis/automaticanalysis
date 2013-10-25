function mask2outline(Mimg)

% Load image first
V = spm_vol(Mimg);
Y = spm_read_vols(V);

% First make an outline of the image
% (requires Image Processing Toolbox)
Y = bwperim(Y);

% Now, let's remove anything that is on the borders of the image
fY = zeros(size(Y));
fY(2:end-1, 2:end-1, 2:end-1) = 1;

% Remove voxels on the border...
Y = Y .* fY;

% Write image back...
spm_write_vol(V,Y);