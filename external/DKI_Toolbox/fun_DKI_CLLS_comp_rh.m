function [D11,D22,D33,D12,D13,D23,...
          W1111, W2222, W3333, W1112, W1113,...
          W1222, W2223, W1333, W2333, W1122,...
          W1133, W2233, W1123, W1223, W1233, S0]=...
fun_DKI_CLLS_comp_rh(data_in,data_mask,bval,bvec)

%                  W1111 W2222 W3333 W1112 W1113
%                  W1222 W2223 W1333 W2333 W1122
%                  W1133 W2233 W1123 W1223 W1233

%% CLLS DKI
% Version implemented by Rafael Neto Henriques for DKIu
% 2/07/2013
%%

Kim=-3/7;
smallc=3;

[Nx, Ny, Nz, Nvol]=size(data_in);

% Minimize ||AX-B||^2
% Where A, B are:

Z15=zeros(Nvol,15);% old Z
Z6=zeros(Nvol,6);% old Z2
Ad=Z6;
Ak=Z15;
Aksbv=Z15;
Adsbv=Z6;

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

% sem b-vals
for v=1:Nvol
    Adsbv(v,1:6)=[bvec(1,v)^2, bvec(2,v)^2, bvec(3,v)^2, ...
        2*bvec(1,v)*bvec(2,v), 2*bvec(1,v)*bvec(3,v), 2*bvec(2,v)*bvec(3,v)];
    
    Aksbv(v,1:15)=[bvec(1,v)^4,...
        bvec(2,v)^4, ...
        bvec(3,v)^4, ...
        4*bvec(1,v)^3*bvec(2,v),...
        4*bvec(1,v)^3*bvec(3,v),...
        4*bvec(2,v)^3*bvec(1,v),...
        4*bvec(2,v)^3*bvec(3,v),...
        4*bvec(3,v)^3*bvec(1,v),...
        4*bvec(3,v)^3*bvec(2,v),...
        6*bvec(1,v)^2*bvec(2,v)^2,...
        6*bvec(1,v)^2*bvec(3,v)^2,...
        6*bvec(2,v)^2*bvec(3,v)^2,...
        12*bvec(1,v)^2*bvec(2,v)*bvec(3,v),...
        12*bvec(2,v)^2*bvec(1,v)*bvec(3,v),...
        12*bvec(3,v)^2*bvec(1,v)*bvec(2,v)];
    
end

%A
A=[ -Ad 1/6*Ak ones(Nvol, 1)]; %SO is now the last column

Z0=zeros(Nx,Ny,Nz);
D11=Z0;
D22=Z0;
D33=Z0;
D12=Z0;
D13=Z0;
D23=Z0;
W1111=Z0;
W2222=Z0;
W3333=Z0;
W1112=Z0;
W1113=Z0;
W1222=Z0;
W2223=Z0;
W1333=Z0;
W2333=Z0;
W1122=Z0;
W1133=Z0;
W2233=Z0;
W1123=Z0;
W1223=Z0;
W1233=Z0;
S0=Z0;

% Ali Tabesh et al fit
bmax=max(bval);
C=zeros(3*Nvol,22);
C(1:Nvol,1:6)=-Adsbv;
C(Nvol+1:2*Nvol, 7:21)=-Aksbv;
C(2*Nvol+1:3*Nvol, 1:6)=-smallc/bmax*Adsbv; C(2*Nvol+1:3*Nvol, 7:21)=Aksbv;

options.MaxIter=100000;
options.LargeScale='off';

d=zeros(3*Nvol,1);

for k=1:Nz
    for j=1:Ny
        for i=1:Nx  
            if(data_mask(i,j,k)==1)  
                %B
                B=log(squeeze(data_in(i,j,k,:)));
                
                indinf=isinf(B);
                if sum(indinf)~=0
                    %minnotinf=min(B(~indinf)); 
                    %maximum diffusion posible take it other 
                    %direction but perhaps other aproaches 
                    %can be better
                    B(indinf)=-10000;  
                end
                
                % ULLS
                piA=pinv(A); %piA pseudoinverse of A
                X0=piA*B;
                
                D0=Adsbv*X0(1:6);
                d(Nvol+1:2*Nvol)=-Kim*D0.*D0;
                
                %% nonlinear fit
                % [X,aa,bb,sol]=lsqcurvefit(@DKImodel,X0,xdata,ydata, lb, ub, options);
                %[X,aa,bb,sol]=lsqcurvefit(@DKImodel,X0,xdata,ydata);
                X=lsqlin(A,B,C,d,[],[],[],[],X0,options);
                
                %% Diffusion parameters
                % diffusion tensor
                D11(i,j,k)=X(1);
                D22(i,j,k)=X(2);
                D33(i,j,k)=X(3);
                D12(i,j,k)=X(4);
                D13(i,j,k)=X(5);
                D23(i,j,k)=X(6);
                
                %mean diffusion
                MD=(D11(i,j,k)+D22(i,j,k)+D33(i,j,k))/3.0;
                
                %% t2 imaging
                S0(i,j,k)=exp(X(22));
                
                %% Kurtosis parameters
                
                W=X(7:21)/(MD^2);
                
                W1111(i,j,k)=W(1);
                W2222(i,j,k)=W(2);
                W3333(i,j,k)=W(3);
                W1112(i,j,k)=W(4);
                W1113(i,j,k)=W(5);
                W1222(i,j,k)=W(6);
                W2223(i,j,k)=W(7);
                W1333(i,j,k)=W(8);
                W2333(i,j,k)=W(9);
                W1122(i,j,k)=W(10);
                W1133(i,j,k)=W(11);
                W2233(i,j,k)=W(12);
                W1123(i,j,k)=W(13);
                W1223(i,j,k)=W(14);
                W1233(i,j,k)=W(15);
            end
        end
    end
    disp(k)
end
