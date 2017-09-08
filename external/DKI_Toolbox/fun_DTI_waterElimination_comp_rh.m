function [D11, D22, D33, D12, D13, D23, S0, F]=fun_DTI_waterElimination_comp_rh(data_in,data_mask,bval,bvec)
%% DTI water elimination model
%% 18/11/2015, Rafael Neto Henriques
optionsF = optimset('Display', 'off');
Diso = 3e-3;
fprecision = 0.01;
ninter = 2;

log_res_sig = -100;

bval = bval(:);

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
F=Z0;

for k=1:Nz
    for j=1:Ny
        for i=1:Nx
            if(data_mask(i,j,k)==1)
                si = double(squeeze(data_in(i,j,k,:)));
                B = log(si);
                indinf=isinf(B);
                if sum(indinf)~=0
                    B(indinf)=log_res_sig;
                end
                % Compute DTI ULLS guess
                X0=A\B;
                
                % Compute DTI WLLS guess
                S = diag(A * X0);
                X0 = pinv(A' * S^2 * A) * A' * S^2 * B;
                
                % WLLS procedure for the water elimination model
                sz = exp(X0(7));
                expbDv = exp(-bval * Diso);
                
                df = 1;
                fin = 0;
                ffn = 1;
                while df > fprecision
                    df = df * 0.1;
                    fw = (fin+df):df:(ffn-df);
                    n = length(fw);
                    expbD = repmat(expbDv, 1, n);
                    [FW, SI] = meshgrid(fw, si);
                    for ii=1:ninter
                        decoupled_si = (SI - (1 - FW)*sz.*expbD)./FW;
                        decoupled_si(decoupled_si<0) = 10^(log_res_sig);
                        Y = log(decoupled_si);
                        S = diag(A * X0);
                        X0 = pinv(A' * S^2 * A) * A' * S^2 * Y;
                        sz = exp(X0(7, :));
                        sz = repmat(sz, Nvol, 1);
                        SIpred = (FW.*exp(A*X0) + (1-FW).*sz.*expbD);
                        Fx2 = sum((SI - SIpred).^2);
                        [Fmin, Fidx] = min(Fx2);
                        f = fw(Fidx);
                        X0 = X0(:, Fidx);
                        sz= exp(X0(7));
                    end
                    fin = f - df;
                    ffn = f + df;
                end
                X = lsqcurvefit(@DTImodel_we, [X0; f], A, si, [0, 0, 0, -3e-3, -3e-3, -3e-3, -100, 0], [3e-3, 3e-3, 3e-3, 3e-3, 3e-3, 3e-3, 100, 1], optionsF);
                
                D11(i,j,k)=X(1);
                D22(i,j,k)=X(2);
                D33(i,j,k)=X(3);
                D12(i,j,k)=X(4);
                D13(i,j,k)=X(5);
                D23(i,j,k)=X(6);
                S0(i,j,k)=exp(X(7));
                F(i,j,k)=X(8);
            end
        end
        disp([k, j]);
    end
end

