function [MK, FA, MD,  l1, l2, l3,v1,v2,v3, Ka, Kr, S0]=fun_DKI_OLS_rh(data_in,data_mask,bval,bvec)

%% ULLS DKI
%
% Neto Henriques R, Correia MM, CamCAN (2012), ‘Towards optimization of diffusion kurtosis imaging to study brain changes with age’, 29th annual meeting of the European Society for Magnetic Resonance in Medicine and Biology, Lisbon, Portugal.
%
% Implemented by Rafael Henriques
% november 2011, MRC-CBU
%%

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
A=[ -Ad 1/6*Ak ones(Nvol, 1)]; %SO is now the last column

Z0=zeros(Nx,Ny,Nz);
l1=Z0;
l2=Z0;
l3=Z0;
v1=zeros(Nx,Ny,Nz,3);
v2=zeros(Nx,Ny,Nz,3);
v3=zeros(Nx,Ny,Nz,3);
MD=Z0;
FA=Z0;
MK=Z0;
S0=Z0;
Ka=Z0;
Kr=Z0;

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
                end
                
                % ULLS
                piA=pinv(A); %piA pseudoinverse of A
                X=piA*B;
                
                %% Diffusion parameters
                % diffusion tensor
                D=[ X(1) X(4) X(5);
                    X(4) X(2) X(6);
                    X(5) X(6) X(3)];
                                          
                [V,L]=eig(D);
                
                %diffusion eigenvalues
                l1(i,j,k)=L(3,3);
                l2(i,j,k)=L(2,2);
                l3(i,j,k)=L(1,1);
                
                %diffusion eigenvectores
                v1(i,j,k,:)=V(:,3);
                v2(i,j,k,:)=V(:,2);
                v3(i,j,k,:)=V(:,1);
                
                %mean diffusion
                MD(i,j,k)=(l1(i,j,k)+l2(i,j,k)+l3(i,j,k))/3.0;
                
                %diffusion fraccional anisotropy
                if(MD(i,j,k)>0)
                    FA(i,j,k)=sqrt(3/2*((l1(i,j,k)-MD(i,j,k))^2+(l2(i,j,k)-MD(i,j,k))^2+(l3(i,j,k)-MD(i,j,k))^2 )/(l1(i,j,k)^2+l2(i,j,k)^2+l3(i,j,k)^2) );
                else
                    FA(i,j,k)=0.0;
                end
                
                % S0
                S0(i,j,k)=exp(X(22));
                
                %% Kurtosis parameters
                % kurtosis "tensor"
                Vd=X(7:21);
                VdMAT=Wcons(Vd);
                % MK
                MK(i,j,k)=GMeanKurt(D,VdMAT,bval,bvec,[],2);
                % directional kurtosis
                
                W=VdMAT/(MD(i,j,k)^2);
                
                p=[V(:,3),V(:,2),V(:,1)];
                
                Wr1111=0;
                for il=1:3
                    for jl=1:3
                        for kl=1:3
                            for ll=1:3
                                Wr1111=Wr1111+W(il,jl,kl,ll)*p(il,1)*p(jl,1)*p(kl,1)*p(ll,1);
                            end
                        end
                    end
                end
                
                k1=(MD(i,j,k)/l1(i,j,k))^2*Wr1111;
                
                Wr2222=0;
                for il=1:3
                    for jl=1:3
                        for kl=1:3
                            for ll=1:3
                                Wr2222=Wr2222+W(il,jl,kl,ll)*p(il,2)*p(jl,2)*p(kl,2)*p(ll,2);
                            end
                        end
                    end
                end
                k2=(MD(i,j,k)/l2(i,j,k))^2*Wr2222;
                
                Wr3333=0;
                for il=1:3
                    for jl=1:3
                        for kl=1:3
                            for ll=1:3
                                Wr3333=Wr3333+W(il,jl,kl,ll)*p(il,2)*p(jl,2)*p(kl,2)*p(ll,2);
                            end
                        end
                    end
                end
                k3=(MD(i,j,k)/l3(i,j,k))^2*Wr3333;
                
                Ka(i,j,k)=k1;
                Kr(i,j,k)=(k2+k3)/2;
                
           end
        end
    end
end
