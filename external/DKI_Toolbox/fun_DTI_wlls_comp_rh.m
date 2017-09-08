function [D11, D22, D33, D12, D13, D23, S0]=fun_DTI_wlls_comp_rh(data_in, data_mask, bval, bvec)
%% DTI water elimination model
%% 18/11/2015, Rafael Neto Henriques

[Nx, Ny, Nz, Nvol]=size(data_in);

% Compute "B-matrix"
Ad=zeros(Nvol,6);
for v=1:Nvol
    b=bval(v);
    Ad(v,1:6)=[b*bvec(1,v)^2, b*bvec(2,v)^2, b*bvec(3,v)^2, ...
        2*b*bvec(1,v)*bvec(2,v), 2*b*bvec(1,v)*bvec(3,v), 2*b*bvec(2,v)*bvec(3,v)];
end
A=[ -Ad ones(Nvol, 1)];

% Initialize variables
Z0=zeros(Nx,Ny,Nz);
D11=Z0;
D22=Z0;
D33=Z0;
D12=Z0;
D13=Z0;
D23=Z0;
S0=Z0;

for k=1:Nz
    for j=1:Ny
        for i=1:Nx
            if(data_mask(i,j,k)==1)
                si = double(squeeze(data_in(i,j,k,:)));
                B = log(si);
                indinf=isinf(B);
                if sum(indinf)~=0
                    B(indinf)=-100;  
                end
                % ULLS
                piA=pinv(A); 
                X0=piA*B;
                
                % Compute DTI WLLS guess
                S = diag(A * X0);
                X = pinv(A' * S^2 * A) * A' * S^2 * B;

                D11(i,j,k)=X(1);
                D22(i,j,k)=X(2);
                D33(i,j,k)=X(3);
                D12(i,j,k)=X(4);
                D13(i,j,k)=X(5);
                D23(i,j,k)=X(6);
                S0(i,j,k)=exp(X(7));
            end
        end
        disp(j)
    end
end
