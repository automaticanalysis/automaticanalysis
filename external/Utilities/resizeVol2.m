function N = resizeVol2(V1,V2)
%%% This function is based on the spm_mask.m function.  V1 is the volume
%%% being resized, V2 is the target Volume being resized to;


% P1 = '/autofs/space/plato_002/users/Aaron/SPM/spm8/apriori/avg152T1_ventricles_MNI.nii';
% P2 = '/autofs/space/plato_002/users/R01_Long_SPM8/RLM055C4_MD_3m/restingState/ss_nn_rr_st_restingState_1.nii';
% V1=spm_vol(P1);
% V2=spm_vol(P2);

% m1=numel(V1);
% m2=numel(V2);

M   = V2(1).mat;
dim = V2(1).dim(1:3);

N = zeros(V2(1).dim);
for j=1:dim(3)
    msk = zeros(dim(1:2));
    Mi  = spm_matrix([0 0 j]);
    M1  = M\V1.mat\Mi;
    img = spm_slice_vol(V1,M1,dim(1:2),[0 NaN]);   
%     keyboard;
%     msk = msk + (img~=0 & isfinite(img));
    msk = msk+img;
    N(:,:,j) = msk;
end

N(N==0) = NaN;

