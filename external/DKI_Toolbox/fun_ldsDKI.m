function fun_ldsDKI(inputDWIfold,inputbvalfold,inputMaskfold)

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

% Computing MK and MD from DKI direct fit
[MK, MD, S0]=fun_DKI_dMK_linear_rh(DWI_data,DWI_mask,bval);

Fit_name=[pth, sc, nnname, '_ldfDKI_MK.nii'];
MKnii=V_mask;
MKnii.hdr.dime.dim(5)=1;
MKnii.hdr.dime.pixdim(5)=0;
MKnii.hdr.dime.xyz_t=0;
MKnii.hdr.dime.datatype=16;
MKnii.img=MK;
save_untouch_nii(MKnii,Fit_name)

Fit_name=[pth, sc, nnname, '_ldfDKI_MD.nii'];
MDnii=V_mask;
MDnii.hdr.dime.dim(5)=1;
MDnii.hdr.dime.pixdim(5)=0;
MDnii.hdr.dime.xyz_t=0;
MDnii.hdr.dime.datatype=16;
MDnii.img=MD;
save_untouch_nii(MDnii,Fit_name)

Fit_name=[pth, sc, nnname, '_ldfDKI_S0.nii'];
S0nii=V_mask;
S0nii.hdr.dime.dim(5)=1;
S0nii.hdr.dime.pixdim(5)=0;
S0nii.hdr.dime.xyz_t=0;
S0nii.hdr.dime.datatype=16;
S0nii.img=MD;
save_untouch_nii(S0nii,Fit_name)
                