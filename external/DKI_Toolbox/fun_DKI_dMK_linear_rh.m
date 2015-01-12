function [MK, MD,  S0]=fun_DKI_dMK_linear_rh(data_in,data_mask,bval)
%% ULLS DKI
% Implemented by Rafael Henriques
% november 2011, MRC-CBU
%%

[Nx, Ny, Nz, Nvol]=size(data_in);


% Minimize ||AX-B||^2
% Where A, B are:

A=[(-bval)' ((bval').^2)/6 ones(Nvol, 1)];

Z0=zeros(Nx,Ny,Nz);
MD=Z0;
MK=Z0;
S0=Z0;
BINF=zeros(Nx,Ny,Nz);

for k=1:Nz
    for j=1:Ny
        for i=1:Nx  
            if(data_mask(i,j,k)==1)  
                %B
                B=log(squeeze(data_in(i,j,k,:)));
                
                indinf=isinf(B);
                if sum(indinf)~=0
                    minnotinf=min(B(~indinf)); 
                    %maximum diffusion posible take it other 
                    %direction but perhaps other aproaches 
                    %can be better
                    B(indinf)=minnotinf;  
                    BINF(i,j,k)=1;
                end
                
                % ULLS
                piA=pinv(A); %piA pseudoinverse of A
                X=piA*B;
                
                %% parameters
                MD(i,j,k)=X(1);         % diffusion
                %V=X(2);                 % Kurtosis
                MK(i,j,k)=X(2)/(MD(i,j,k)^2);
                S0(i,j,k)=exp(X(3));    % t2 imaging
     
            end
        end
    end
end
