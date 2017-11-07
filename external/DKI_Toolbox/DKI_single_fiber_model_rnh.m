function [AWF,Da,De,tortu,ADa,ADe,RDe,Kmax]=...
    DKI_single_fiber_model_rnh(D,K,V,Nvis,Uvis,optionsT)
% This function computes biophysical measures from DKI based on the
% assumption that the are working with white matter fibers that are well
% aligned and represented by two media - the intra and extra cellular media
%
% Inputs:
% Dv=[D11 D22 D33 D12 D13 D23];
% Kv=[W1111 W2222 W3333 W1112 W1113...
%     W1222 W2223 W1333 W2333 W1122...
%     W1133 W2233 W1123 W1223 W1233];
% V = Initial search grid
% Nvis = number of grid point neighbour 
% Uvis = neighbour index
% optionsT = convergence options
% 
% Output:
% AWF = axonal water fraction
% Da = diffusion tensor of the intracelular component model 
% De = diffusion tensor of the extracelular component model
% tortu = tortuosity
%
% Metrics suggested by Fieremans et al. (2011),
% White matter characterization with diffusional kurtosis imaging
% Implemented by Rafael Neto Henriques 2013
% Last update 11/04/2015
% New convergence procedures

% Step 1 sampling initial grid
RK=DirectionalKurt_V(D,K,V);

% step 2 - find maxima of KT
np=length(RK);
MeM=false(np,1);
for vt=1:np
    condi=(RK(vt)>RK(Uvis(vt,1:Nvis(vt))));
    if sum(condi)==Nvis(vt)
        MeM(vt)=true;
    end
end

% step 3 - select the maxima points
MR=RK(MeM);
MV=V(MeM,:);

% step 4 - compute metrics case of spherical KT
if isempty(MR) 
    AWF=0;
    Da=zeros(1,6);
    De=Da;
    tortu=0;
    ADa=0;
    ADe=0;
    RDe=0;
    Kmax=0;
else
    % step 5 - convergence of maxima for non-spherical KT
    MD=mean(D(1:3));
    [az,el] = cart2sph(MV(:,1),MV(:,2),MV(:,3));
    lMR=length(MR);
    for vm=1:lMR
        xi(2)=el(vm);
        xi(1)=az(vm);
        myf=@(x) DirectionalKurt(D,K,x,true,MD);
        xf = fminsearch(myf,xi,optionsT);
        az(vm)=xf(1);
        el(vm)=xf(2);
        MR(vm)=DirectionalKurt(D,K,xf,false,MD);
    end
    Kmax=max(MR);
    
    % step 6 - compute AWF 
    if Kmax > 12
        Kmax=12;
        AWF=0.8;
    else
        AWF=Kmax/(Kmax+3);
    end
    % step 7 - compute intra and extracelular tensors Da and De
    [De,Da,DTe,DTa]=EQ10EQ11_Fieremans(D,K,AWF);
    
    % step 8 - compute Da and De metrics 
    [ADa,ADe,RDe,tortu]=EQ12EQ13EQ14EQ15_Fieremans(DTe,DTa);
end