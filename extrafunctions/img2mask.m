function img2mask(Mimg)

% Load image first
V = spm_vol(Mimg);
Y = spm_read_vols(V);

% First round the image...
Y = round(Y);

% Anything above 1 is 1
Y(Y>1) = 1;

% Any NaNs are 0
Y(isnan(Y)) = 0;

% Write image back...
spm_write_vol(V,Y);