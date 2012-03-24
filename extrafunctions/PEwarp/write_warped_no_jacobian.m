function write_warped_no_jacobian(x, xscale, order, source, warped)

sourceVol = spm_vol(source);

[gridX, gridY, gridZ] = ndgrid(1 : sourceVol.dim(1), 1 : sourceVol.dim(2), 1 : sourceVol.dim(3));
bx = spm_dctmtx(sourceVol.dim(1), order(1));
by = spm_dctmtx(sourceVol.dim(2), order(2));
bz = spm_dctmtx(sourceVol.dim(3), order(3));

def = spm_get_def(bx, by, bz, x * xscale);
warpedData = reshape(spm_sample_vol(sourceVol, gridX(:), gridY(:) + def, gridZ(:), 1), sourceVol.dim);

warpedVol = sourceVol;
warpedVol.fname = warped;

spm_write_vol(warpedVol, warpedData);
