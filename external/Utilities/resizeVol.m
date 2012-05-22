function N = resizeVol(V1,V2)
%%% This function is based on the spm_mask.m function.  V1 is the volume
%%% being resized, V2 is the target Volume being resized to;
%%% 
%%% Note: The inputs are the structures returned from spm_vols.
%%%
%%% Written by Aaron Schultz, May 5th, 2010 (aschultz@martinos.org)

M   = V2(1).mat;
dim = V2(1).dim(1:3);

N = zeros(V2(1).dim);
for j=1:dim(3)
    msk = zeros(dim(1:2));
    Mi  = spm_matrix([0 0 j]);
    M1  = M\V1.mat\Mi;
    img = spm_slice_vol(V1,M1,dim(1:2),[0 NaN]);    
    msk = msk + (img~=0 & isfinite(img));
%     msk = msk+img;
    N(:,:,j) = msk;
end

N(N==0) = NaN;

