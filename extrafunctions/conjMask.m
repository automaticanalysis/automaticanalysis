%% conjMask creates a mask image that includes voxels common to both
% Importantly, the mask can only work on images of same size
% This function places the resulting image in the directory of Img1
% and names it similarly to Img1
% Also, you may add a prefix of your choice to the final image
function conjMask(Img1, Img2, thresh, prefix)

if nargin < 3 || isempty(thresh)
    thresh = [0 0];
end
if nargin < 4 || isempty(prefix)
    prefix = 'c';
    aas_log([],false,'Using default prefix ''c''')
end

V1 = spm_vol(Img1);
V2 = spm_vol(Img2);

Y1 = spm_read_vols(V1);
Y2 = spm_read_vols(V2);

Vc = V1;
[pth nme ext] = fileparts(Vc.fname);
Vc.fname = fullfile(pth, [prefix nme ext]);
try
Yc = double(and((Y1>thresh(1)),(Y2>thresh(2))));

spm_write_vol(Vc, Yc);
if ~(sum(Yc(:)>0) > 0)
    aas_log([],false,'Failed not make conjunction: no overlap')
end

[junk, nImg1] = fileparts(Img1);
[junk, nImg2] = fileparts(Img2);
aas_log([],false,sprintf('Voxel count: \n\t%s: %.0f\n\t%s: %.0f\n\tConjunction: %.0f', ...
    nImg1, sum(Y1(:)>thresh(1)), nImg2, sum(Y2(:)>thresh(2)), sum(Yc(:)>0)))
catch
    if ~all(size(Y1) == size(Y2));
        aas_log([],true,'Failed to make conjunction: different sized files')
    else
        aas_log([],true,'Failed to make conjunction: UNKNOWN REASON!')
    end
end