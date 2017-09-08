path = 'C:\Users\Rafael\Desktop\Subject1\';

addpath('C:\Users\Rafael\Desktop\DKIu_v1.2\NIFTI_toolbox\')

z = 30;

figure
DataDWI = [path, 'DTImetrics/DWI_brain_GF150_olsDKI_DT_MD.nii'];
V = load_nii(DataDWI, [], [], [], [], [], 1);
MK = double(V.img);
trans = squeeze(MK(:,:,z));
trans = trans';
transc = trans;
trans = transc(end:-1:1,:);
imagesc(trans, [0 2e-3])
axis off
colormap gray

figure
DataDWI = [path, 'DTImetrics/DWI_brain_GF150_olsDKI_DT_FA.nii'];
V = load_nii(DataDWI, [], [], [], [], [], 1);
MK = double(V.img);
trans = squeeze(MK(:,:,z));
trans = trans';
transc = trans;
trans = transc(end:-1:1,:);
imagesc(trans, [0 1.5])
axis off
colormap gray

figure
DataDWI = [path, 'DKImetrics/DWI_brain_GF150_olsDKI_KT_MK.nii'];
V = load_nii(DataDWI, [], [], [], [], [], 1);
MK = double(V.img);
trans = squeeze(MK(:,:,z));
trans = trans';
transc = trans;
trans = transc(end:-1:1,:);
imagesc(trans, [0 1.5])
axis off
colormap gray