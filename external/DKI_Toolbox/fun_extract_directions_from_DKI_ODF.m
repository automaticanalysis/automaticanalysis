function [NM, MR, FDir] = fun_extract_directions_from_DKI_ODF(DT,KT,Mask,alfa)
% Rafael Neto Henriques
% Implemented 23/04/2014
% Updates 01/09/2014 - new convergence procedure based on Matlab's
% quasi-Newton algorithm
% Updates 07/04/2015 - New and faster subfunctions (4 times faster)

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
end