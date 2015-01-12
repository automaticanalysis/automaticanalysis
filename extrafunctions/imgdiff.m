function imgdiff(Mimg, Simg, Dimg)

if nargin < 3
    [Mpth Mfn Mext] = fileparts(Mimg);
    Dimg = fullfile(Mpth, ['diff_' Mfn Mext]);
end

% Load images first
mV = spm_vol(Mimg);
mY = spm_read_vols(mV);

sV = spm_vol(Simg);
sY = spm_read_vols(sV);

dV = mV;
dY = mY - sY;
dV.fname = Dimg;
spm_write_vol(dV, dY);
