function [DT, KT, S0]=fun_DKI_RLS_rh(data_in,data_mask,bval,bvec,DTi,KTi,S0i,RStd)
% DKI rician fit
% Implemented by Rafael Neto Henriques, MRC-CBU April 2014
%
% This riqueres a first estimate of DT abd KT tensor
% 
% DT and KT tensors elements have to be in this order:
% DT elements:     D11 D22 D33 D12 D13 D23
% KT elements:
%                  W1111 W2222 W3333 W1112 W1113
%                  W1222 W2223 W1333 W2333 W1122
%                  W1133 W2233 W1123 W1223 W1233

[Nx, Ny, Nz, Nvol]=size(data_in);

% Minimize ||AX-B||^2
% Where A, B are:

Z15=zeros(Nvol,15);% old Z
Z6=zeros(Nvol,6);% old Z2
Ad=Z6;
Ak=Z15;

for v=1:Nvol
    b=bval(v);
    Ad(v,1:6)=[b*bvec(1,v)^2, b*bvec(2,v)^2, b*bvec(3,v)^2, ...
        2*b*bvec(1,v)*bvec(2,v), 2*b*bvec(1,v)*bvec(3,v), 2*b*bvec(2,v)*bvec(3,v)];
    
    Ak(v,1:15)=[b*b*bvec(1,v)^4,...
        b*b*bvec(2,v)^4, ...
        b*b*bvec(3,v)^4, ...
        4*b*b*bvec(1,v)^3*bvec(2,v),...
        4*b*b*bvec(1,v)^3*bvec(3,v),...
        4*b*b*bvec(2,v)^3*bvec(1,v),...
        4*b*b*bvec(2,v)^3*bvec(3,v),...
        4*b*b*bvec(3,v)^3*bvec(1,v),...
        4*b*b*bvec(3,v)^3*bvec(2,v),...
        6*b*b*bvec(1,v)^2*bvec(2,v)^2,...
        6*b*b*bvec(1,v)^2*bvec(3,v)^2,...
        6*b*b*bvec(2,v)^2*bvec(3,v)^2,...
        12*b*b*bvec(1,v)^2*bvec(2,v)*bvec(3,v),...
        12*b*b*bvec(2,v)^2*bvec(1,v)*bvec(3,v),...
        12*b*b*bvec(3,v)^2*bvec(1,v)*bvec(2,v)];
    
end

%A
A=[ -Ad/1000 1/6*Ak/1000000 ones(Nvol, 1)]; %SO is now the last column

DT=zeros(Nx, Ny, Nz, 6);
KT=zeros(Nx, Ny, Nz, 15);
S0=zeros(Nx, Ny, Nz);

for k=1:Nz
    for j=1:Ny
        for i=1:Nx  
            if(data_mask(i,j,k)==1)  
                y=squeeze(data_in(i,j,k,:));
                DT0kji=squeeze(DTi(i,j,k,:));
                KT0kji=squeeze(KTi(i,j,k,:));
                S0kji=squeeze(S0i(i,j,k));
                [DTkji,KTkji,S0kji]=DKI_likelihood_rnh(DT0kji,KT0kji,S0kji,A,y,RStd);
                DT(i,j,k,:)=DTkji;
                KT(i,j,k,:)=KTkji;
                S0(i,j,k)=S0kji;
                disp(i)
            end
        end
        disp(j)
    end
    disp(k)
end
