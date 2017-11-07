clear all
close all
optionsF = optimset('Display', 'off');
Diso = 3e-3;

% White matter component
DAi = 1.7e-3;
DRi = 0.3e-3;
fie = 0.7;
fp = 1;

% Water component
DAe = Diso;
DRe = Diso;

% Other parameters
fphi = [0 0];
ftheta = [0 0];
S0=200;
sig=0;

% Scanner parameters (2 b-values, 125 directions and 2 b0)
load Dir125.mat
X = V;
bval=[500*ones(125,1); 750*ones(125,1); 1000*ones(125,1); 1500*ones(125,1); 0; 0];
bvec=[X; X; X; X; 0,0,0; 0,0,0]';

% Run simulations
si=Sgenerator(DAi,DRi,DAe,DRe,fie,fp,fphi,ftheta,bvec,bval,S0,sig);

% algorithms parameters
fprecision = 0.01;

% Compute B-matrix
Nvol=length(si);
Z6=zeros(Nvol,6);% old Z2
Ad=Z6;
for v=1:Nvol
    b=bval(v);
    Ad(v,1:6)=[b*bvec(1,v)^2, b*bvec(2,v)^2, b*bvec(3,v)^2, ...
        2*b*bvec(1,v)*bvec(2,v), 2*b*bvec(1,v)*bvec(3,v), 2*b*bvec(2,v)*bvec(3,v)];
end
A=[ -Ad ones(Nvol, 1)];

tic
% Compute DTI first guess
B=log(si);
piA=pinv(A);
X0=A\B;
% Compute DTI WLLS guess
S = diag(A * X0);
X0 = pinv(A' * S^2 * A) * A' * S^2 * B;
Xdti = X0;
%WLLS procedure for the water elimination model
S0 = exp(X0(7));
Diso = 3e-3;
expbDv = exp(-bval * Diso);

df = 1;
fin = 0;
ffn = 1;
ninter = 2;
while df > fprecision
    df = df * 0.1;
    fw = (fin+df):df:(ffn-df);
    n = length(fw);
    expbD = repmat(expbDv, 1, n);
    [FW, SI] = meshgrid(fw, si);
    for i=1:ninter
        Y = log((SI - (1 - FW)*S0.*expbD)./FW);
        S = diag(A * X0);
        X0 = pinv(A' * S^2 * A) * A' * S^2 * Y;
        S0 = exp(X0(7, :));
        S0 = repmat(S0, Nvol, 1);
        
        SIpred = (FW.*exp(A*X0) + (1-FW).* S0.*expbD);
        Fx2 = sum((SI - SIpred).^2);
        [Fmin, Fidx] = min(Fx2);
        f = fw(Fidx);
        X0 = X0(:, Fidx);
        S0 = exp(X0(7));
    end
    fin = f - df;
    ffn = f + df;
end
disp(toc)
X0(1:6)=X0(1:6)*1000;
X0 = [X0; f];
A(:, 1:6) = A(:, 1:6)/1000;
X = lsqcurvefit(@DTImodel_we, X0, A, si, [], [], optionsF);
X(1:6) = X(1:6)/1000;
f = X(8);
disp(toc)