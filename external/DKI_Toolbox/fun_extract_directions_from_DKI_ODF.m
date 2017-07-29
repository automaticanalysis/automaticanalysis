function fun_extract_directions_from_DKI_ODF(DT_name,KT_name,Mask_name,alfa)
% Rafael Neto Henriques
% Implemented 23/04/2014
% Updates 01/09/2014 - new convergence procedure based on Matlab's
% quasi-Newton algorithm
% Updates 07/04/2015 - New and faster subfunctions (4 times faster)

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
disp('Processing maxima DKI-ODF directions ...')
for k=1:N3
    for j=1:N2
        for i=1:N1
            if(Mask(i,j,k)==1)
                D=DT(i,j,k,:);
                
                W=KT(i,j,k,:);
                % Findin directions from the DKI based ODF
                [MVi,MRi]=Procedure_DKIODF(D,W,V,Nvis,Uvis,alfa,optionsT);
                
                % DKI-ODF is know to resolve three fiber direction.
                % Therefore when algorithm gives more than three fiber we
                % revome the suspious small amplitude fibers
                NMi=length(MRi);
                if NMi>3
                    [sMRi,iMRi]=sort(MRi,'descend');
                    MRi=sMRi(1:3);
                    MVi=MVi(iMRi(1:3),:);
                    NMi=3;
                    NM(i,j,k)=NMi;
                    MR(i,j,k,1:NMi)=MRi;
                    FDir(i,j,k,1:NMi,:)=MVi;
                elseif NMi>0
                    [sMRi,iMRi]=sort(MRi,'descend');
                    MRi=sMRi;
                    MVi=MVi(iMRi,:);
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
output=[pth,'/DKI_tract_m/', outname, '_DKIODF_alfa',num2str(alfa),'_NM.nii'];
Vol=V_Mask;
Vol.hdr.dime.datatype=16;
Vol.img=NM;
save_untouch_nii(Vol,output);

%MR
output=[pth,'/DKI_tract_m/', outname, '_DKIODF_alfa',num2str(alfa),'_MR.nii'];
Vol.hdr.dime.dim(5)=3;
Vol.img=MR;
save_untouch_nii(Vol,output);

%FDir
output=[pth,'/DKI_tract_m/', outname, '_DKIODF_alfa',num2str(alfa),'_FDir.nii'];
Vol.hdr.dime.dim(6)=3;
Vol.img=FDir;
save_untouch_nii(Vol,output);

end