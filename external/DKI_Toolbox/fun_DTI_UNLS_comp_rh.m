function [D11,D22,D33,D12,D13,D23, S0]=fun_DTI_UNLS_comp_rh(data_in,data_mask,bval,bvec)
%% UNLS DTI
% First version implemented by Rafael Henriques
% november 2011, MRC-CBU
% comp version 4 July 2013
%%
optionsF = optimset('Display', 'off');

[Nx, Ny, Nz, Nvol]=size(data_in);

% Minimize ||AX-B||^2
% Where A, B are:

Z6=zeros(Nvol,6);% old Z2
Ad=Z6;


for v=1:Nvol
    b=bval(v);
    Ad(v,1:6)=[b*bvec(1,v)^2, b*bvec(2,v)^2, b*bvec(3,v)^2, ...
        2*b*bvec(1,v)*bvec(2,v), 2*b*bvec(1,v)*bvec(3,v), 2*b*bvec(2,v)*bvec(3,v)];  
end

%A
A=[ -Ad ones(Nvol, 1)]; %SO is now the last column
xdata=A';

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
                %B
                Sig=double(squeeze(data_in(i,j,k,:)));
                indzero=(Sig==0);
                if sum(indzero)>0
                    Sig(indzero)=eps;
                end
                B=log(Sig);
                ydata=Sig'; % problema dos B ja nao e' totalmente necessario
 
                
                %% First gest ULLS
                piA=pinv(A); %piA pseudoinverse of A
                X0=piA*B;
                
                X0(1:6)=X0(1:6)*1000;
                
                %% nonlinear fit
                % [X,aa,bb,sol]=lsqcurvefit(@DKImodel,X0,xdata,ydata, lb, ub, options);
                %[X,aa,bb,sol]=lsqcurvefit(@DKImodel,X0,xdata,ydata);
                X=lsqcurvefit(@DTImodel2,X0,xdata,ydata,[],[],optionsF);
                
                X(1:6)=X(1:6)/1000;
                
                %% Diffusion parameters
                % diffusion tensor
                
                D11(i,j,k)=X(1);
                D22(i,j,k)=X(2);
                D33(i,j,k)=X(3);
                D12(i,j,k)=X(4);
                D13(i,j,k)=X(5);
                D23(i,j,k)=X(6);
                S0(i,j,k)=exp(X(7));
                
            end
        end
    end
    fprintf('%d of %d slices done!\n',k,Nz);
end
