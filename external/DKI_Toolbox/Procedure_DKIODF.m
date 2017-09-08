function [MV,MR,R]=Procedure_DKIODF(Dv,Kv,V,Nvis,Uvis,alfa,optionsT)
% Implemented on 2/07/2014
% Update 6/04/2015

% Inputs:
% Dv=[D11 D22 D33 D12 D13 D23];
% Kv=[W1111 W2222 W3333 W1112 W1113...
%     W1222 W2223 W1333 W2333 W1122...
%     W1133 W2233 W1123 W1223 W1233];
% First
% Output:
% MV = vertices correspondente a maximos
% MR = valores dos vertices maximos
% R = valores dos vertices todos

% step 0 - difine parametres to speed DKIODF computation
MD=mean(Dv(1:3));
DT=[Dv(1),Dv(4),Dv(5);...
    Dv(4),Dv(2),Dv(6);...
    Dv(5),Dv(6),Dv(3) ];
U=pinv(DT)*MD;


% step 1 - Find the values for the initial points
R=DKIODF_v(Dv,Kv,V,alfa);

% step 2 - detect maxima within the initial points
np=length(R);
MeM=false(np,1);
for vt=1:np
    condi=(R(vt)>R(Uvis(vt,1:Nvis(vt))));
    if sum(condi)==Nvis(vt)
        MeM(vt)=true;
    end
end

% step 3 - select the maxima points
MV=V(MeM,:);
MR=R(MeM,:);

% step 4 - convergence
if ~isempty(MR)
    [az,el] = cart2sph(MV(:,1),MV(:,2),MV(:,3));
    lMR=length(MR);
    for vm=1:lMR
        xi(2)=el(vm);
        xi(1)=az(vm);
        neg=true;
        myf=@(x) DKIODF(Dv,Kv,x,alfa,neg,U);
        xf = fminsearch(myf,xi,optionsT);
        az(vm)=xf(1);
        el(vm)=xf(2);
        MR(vm)=DKIODF(Dv,Kv,xf,alfa,false,U);
    end
    
    [x,y,z] = sph2cart(az,el,1);
    MV=[x,y,z];
    
    % step 5 select only unique directions
    vaf=(sum(abs(triu(MV*MV',1))>0.99)==0);
    MV=MV(vaf,:);
    MR=MR(vaf);
    [MR,id]=sort(MR,'descend');
    MV=MV(id,:);
    
end
end