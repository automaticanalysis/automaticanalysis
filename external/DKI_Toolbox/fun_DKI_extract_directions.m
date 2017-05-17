function [NM, MR, FDir] = fun_DKI_extract_directions(DT,DK,Mask,alfa)
% Rafael Neto Henriques
% Implemented 23/04/2014
% Last update 01/09/2014 - new convergence procedure based on Matlab's
% quasi-Newton algorithm
% RK based on pK

% Grid search
V = Dir125('V');
Uvis = Dir125('Uvis');
Nvis = Dir125('Nvis');

% Inicialize
[N1,N2,N3]=size(Mask);
NM=zeros(N1,N2,N3); % Number of maximums detected
%RK=zeros(N1,N2,N3); % Number of maximums detected
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

% Convergence parameters
optionsT = optimset('TolX',1e-3);

% Compute orientational parameters
for k=1:N3
    for j=1:N2
        for i=1:N1
            if(Mask(i,j,k)==1)
                D=DT(i,j,k,:);
                
                W=DK(i,j,k,:);
                % Maximos do NG-ODF da lazar (MV) e respectivas orientacoes (MV)
                [MVi,MRi]=fun_DKI_ODFprocess(D,W,V,Nvis,Uvis,alfa,optionsT);
                
                % For my previous analysis the tensor can have local
                % minimuns which does not correspond to fiber directions.
                % Therefore three directions with larger values of kurtosis
                % (maximum directions detected by DKI) are only saved.
                % Lower values should not correspond to fiber directions
                NMi=length(MRi);
                if NMi>3
                    [sMRi,iMRi]=sort(MRi,'descend');
                    MRi=sMRi(1:3);
                    MVi=MVi(iMRi(1:3),:);
                    NMi=3;
                else
                    [sMRi,iMRi]=sort(MRi,'descend');
                    MRi=sMRi;
                    MVi=MVi(iMRi,:);
                end
                
                %                 RKi=0;
                %                 for f=1:NMi
                %                     RKi=RKi+pKT_cart(D,W,MVi(f,1),MVi(f,2),MVi(f,3));
                %                 end
                %                 RKi=RKi/NMi;
                %
                %                 RK(i,j,k)=RKi;
                NM(i,j,k)=NMi;
                MR(i,j,k,1:NMi)=MRi;
                FDir(i,j,k,1:NMi,:)=MVi;
            end
        end
    end
end

end