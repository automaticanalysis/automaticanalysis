function fun_extract_directions_from_DKI_KT(DT_name,KT_name,Mask_name, fMD)
% Rafael Neto Henriques
%
% Updates: 01/09/2014
% New convergence procedure (based on quasi-newton algorithm provided by Matlab)
% Correction based on a setted MD percentage threshold
% Correction not done if fibre direction estimate is only one
% 09/04/2015
% New and faster subfunctions
% RK is the weighted avarege (consistent to article Neto Henriques et al.,
% 2015)

% load data
V_DT=load_untouch_nii(DT_name);
V_KT=load_untouch_nii(KT_name);
V_Mask=load_untouch_nii(Mask_name);

DT=V_DT.img;
KT=V_KT.img;
Mask=V_Mask.img;

% Inicialize

[N1,N2,N3]=size(Mask);
NM=zeros(N1,N2,N3); % Number of maximums detected
RK=zeros(N1,N2,N3); % Number of maximums detected
FDir=zeros(N1,N2,N3,3,3); % Acording to my previous analysis kurtosis 
                          % tensor is limited to distiguish 3 non-complanar
                          % fibers therefore we will extract at least 3
                          % directions per voxel
                          %
                          % forth dimension correspond to different
                          % direntions of fibers (maximum 3)
                          %
                          % fift coordinates x, y, z of each fiber direction
MR=zeros(N1,N2,N3,3); % Kurtosis values for each direction detected and saved

% Grid search
load('Dir125.mat')

% Convergence parameters
optionsT = optimset('TolX',1e-3,'Display', 'off');

% Compute orientational parameters
disp('Processing maxima perpendicular kurtosis directions ... ')
for k=1:N3
    for j=1:N2
        for i=1:N1
            if(Mask(i,j,k)==1)
                D=DT(i,j,k,:);
                
                W=KT(i,j,k,:);
                [MVi,MRi]=Procedure_KT(D,W,V,Nvis,Uvis,optionsT,fMD);
                NMi=length(MRi);
                if NMi>3
                    [sMRi,iMRi]=sort(MRi,'descend');
                    MRi=sMRi(1:3);
                    MVi=MVi(iMRi(1:3),:);
                    NMi=3;
                    RK(i,j,k)=MRi'/sum(MRi)*MRi;
                    NM(i,j,k)=NMi;
                    MR(i,j,k,1:NMi)=MRi;
                    FDir(i,j,k,1:NMi,:)=MVi;
                elseif NMi>0
                    [sMRi,iMRi]=sort(MRi,'descend');
                    MRi=sMRi;
                    MVi=MVi(iMRi,:);
                    RK(i,j,k)=MRi'/sum(MRi)*MRi;
                    NM(i,j,k)=NMi;
                    MR(i,j,k,1:NMi)=MRi;
                    FDir(i,j,k,1:NMi,:)=MVi;
                end
            end
        end
    end
    disp(k)
end

% save
[pth,outname, daa]=fileparts(KT_name);
mkdir([pth,'/DKI_tract_m/']);

%NM
output=[pth,'/DKI_tract_m/', outname, '_DKIKT_MDC',num2str(fMD*100),'_NM.nii'];
Vol=V_Mask;
Vol.hdr.dime.datatype=16;
Vol.img=NM;
save_untouch_nii(Vol,output);

output=[pth,'/DKI_tract_m/', outname, '_DKIKT_MDC',num2str(fMD*100),'_newRK.nii'];
Vol.img=RK;
save_untouch_nii(Vol,output);

%MR
output=[pth,'/DKI_tract_m/', outname, '_DKIKT_MDC',num2str(fMD*100),'_MR.nii'];
Vol.hdr.dime.dim(5)=3;
Vol.img=MR;
save_untouch_nii(Vol,output);

%FDir
output=[pth,'/DKI_tract_m/', outname, '_DKIKT_MDC',num2str(fMD*100),'_FDir.nii'];
Vol.hdr.dime.dim(6)=3;
Vol.img=FDir;
save_untouch_nii(Vol,output);

end