function fun_waterelimination_DTI(inputDWIfold,inputbvalfold,inputbvecfold,inputMaskfold)

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

% non linear DTI fit
[D11,D22,D33,D12,D13,D23, S0, F]=...
    fun_DTI_waterElimination_comp_rh(DWI_data,DWI_mask,bval,bvec);

% save DTI tensor
DT(:,:,:,6)=D23;
DT(:,:,:,1)=D11;
DT(:,:,:,2)=D22;
DT(:,:,:,3)=D33;
DT(:,:,:,4)=D12;
DT(:,:,:,5)=D13;

% For other software compatibility
DTI(:,:,:,6)=D33;
DTI(:,:,:,1)=D11;
DTI(:,:,:,2)=D12;
DTI(:,:,:,3)=D22;
DTI(:,:,:,4)=D13;
DTI(:,:,:,5)=D23;

Fit_name=[pth, sc, nnname, '_waterEliminationDTI_DT.nii'];
DTnii=V_DWI;
DTnii.hdr.dime.dim(5)=6;
DTnii.hdr.dime.pixdim(5)=0;
DTnii.hdr.dime.xyz_t=0;
DTnii.hdr.dime.datatype=16;
DTnii.img=DT;
save_untouch_nii(DTnii,Fit_name)

Fit_name=[pth, sc, nnname, '_waterEliminationDTI_S0.nii'];
S0nii=V_mask;
S0nii.hdr.dime.dim(5)=1;
S0nii.hdr.dime.pixdim(5)=0;
S0nii.hdr.dime.xyz_t=0;
S0nii.hdr.dime.datatype=16;
S0nii.img=S0;
save_untouch_nii(S0nii,Fit_name)

Fit_name=[pth, sc, nnname, '_waterEliminationDTI_F.nii'];
Fnii=V_mask;
Fnii.hdr.dime.dim(5)=1;
Fnii.hdr.dime.pixdim(5)=0;
Fnii.hdr.dime.xyz_t=0;
Fnii.hdr.dime.datatype=16;
Fnii.img=F;
save_untouch_nii(Fnii,Fit_name)

DTInii=V_DWI;
DTInii.hdr.dime.dim(5)=6;
DTInii.hdr.dime.datatype=16;
DTInii.hdr.dime.pixdim(5)=0;
DTInii.hdr.dime.xyz_t=0;
DTInii.img=DTI;
save_untouch_nii(DTInii,[pth, sc, nnname, '_waterEliminationDTI_DTorient.nii'])