function fun_olsDKI(inputDWIfold,inputbvalfold,inputbvecfold,inputMaskfold)

archstr = computer('arch');
st_comp=archstr(1:3);
if strcmp(st_comp,'win')
    sc='\';
else
    sc='/';
end

[pth, nnname, daa]=fileparts(inputDWIfold);

V_DWI=load_untouch_nii(inputDWIfold); DWI_data=V_DWI.img;
V_mask=load_untouch_nii(inputMaskfold); DWI_mask=V_mask.img;
bval=load(inputbvalfold);
bvec=load(inputbvecfold);

% Computing DT and DKI Tensor
[D11,D22,D33,D12,D13,D23,...
    W1111, W2222, W3333, W1112, W1113,...
    W1222, W2223, W1333, W2333, W1122,...
    W1133, W2233, W1123, W1223, W1233, S0]=...
    fun_DKI_ULLS_comp_rh(double(DWI_data),DWI_mask,bval,bvec);

% 22 variaveis do DKI
DT(:,:,:,6)=D23;
DT(:,:,:,1)=D11;
DT(:,:,:,2)=D22;
DT(:,:,:,3)=D33;
DT(:,:,:,4)=D12;
DT(:,:,:,5)=D13;

KT(:,:,:,15)=W1233;
KT(:,:,:,1)=W1111;
KT(:,:,:,2)=W2222;
KT(:,:,:,3)=W3333;
KT(:,:,:,4)=W1112;
KT(:,:,:,5)=W1113;
KT(:,:,:,6)=W1222;
KT(:,:,:,7)=W2223;
KT(:,:,:,8)=W1333;
KT(:,:,:,9)=W2333;
KT(:,:,:,10)=W1122;
KT(:,:,:,11)=W1133;
KT(:,:,:,12)=W2233;
KT(:,:,:,13)=W1123;
KT(:,:,:,14)=W1223;

% For other software compatibility 
DTI(:,:,:,6)=D33;
DTI(:,:,:,1)=D11;
DTI(:,:,:,2)=D12;
DTI(:,:,:,3)=D22;
DTI(:,:,:,4)=D13;
DTI(:,:,:,5)=D23;

Fit_name=[pth, sc, nnname, '_olsDKI_KT.nii'];

KTnii=V_DWI;
KTnii.hdr.dime.dim(5)=15;
KTnii.hdr.dime.pixdim(5)=0;
KTnii.hdr.dime.xyz_t=0;
KTnii.hdr.dime.datatype=16;
KTnii.img=KT;
save_untouch_nii(KTnii,Fit_name)

Fit_name=[pth, sc, nnname, '_olsDKI_DT.nii'];
DTnii=V_DWI;
DTnii.hdr.dime.dim(5)=6;
DTnii.hdr.dime.pixdim(5)=0;
DTnii.hdr.dime.xyz_t=0;
DTnii.hdr.dime.datatype=16;
DTnii.img=DT;
save_untouch_nii(DTnii,Fit_name)

Fit_name=[pth, sc, nnname, '_olsDKI_S0.nii'];
S0nii=V_mask;
S0nii.hdr.dime.dim(5)=1;
S0nii.hdr.dime.pixdim(5)=0;
S0nii.hdr.dime.xyz_t=0;
S0nii.hdr.dime.datatype=16;
S0nii.img=S0;
save_untouch_nii(S0nii,Fit_name)

DTInii=V_DWI;
DTInii.hdr.dime.dim(5)=6;
DTInii.hdr.dime.datatype=16;
DTInii.hdr.dime.pixdim(5)=0;
DTInii.hdr.dime.xyz_t=0;
DTInii.img=DTI;
save_untouch_nii(DTInii,[pth, sc, nnname, '_olsDKI_DTorient.nii'])