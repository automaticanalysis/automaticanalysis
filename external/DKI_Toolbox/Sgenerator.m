function S=Sgenerator(DAi,DRi,DAe,DRe,fie,fp,fphi,ftheta,bvec,bval,S0,sig)

NoF = length(fp);

Di=[DAi 0 0;0 DRi 0;0 0 DRi];
De=[DAe 0 0;0 DRe 0;0 0 DRe];

AllDi=zeros(3,3,NoF);
AllDe=AllDi;
f_p1=fp; % compartment fractions part 1 (intracelular compartments)
f_p2=fp; % compartment fractions part 2 (extracelular compartments)

for c=1:NoF
    [xa,ya,za]=sph2cart(fphi(c),ftheta(c),1);
    Xa=[xa,ya,za];
    
    if fphi(c) == pi/2 && ftheta(c) ==0
        Xb=cross(Xa,[1 0 0]);
        Xb=Xb/norm(Xb);
        Xc=cross(Xa,Xb);
        Xc=Xc/norm(Xc);
    else
        Xb=cross(Xa,[0 1 0]);
        Xb=Xb/norm(Xb);
        Xc=cross(Xa,Xb);
        Xc=Xc/norm(Xc);
    end
    
    f_p1(c)=fp(c)*fie;
    f_p2(c)=fp(c)*(1-fie);
    
    Evector=[Xa',Xb',Xc'];
    
    Dci=Evector*Di*pinv(Evector);
    AllDi(:,:,c)=Dci;
    
    Dce=Evector*De*pinv(Evector);
    AllDe(:,:,c)=Dce;

end

f=[f_p1 f_p2];
AllD=zeros(3,3,NoF*2);
AllD(:,:,1:NoF)=AllDi;
AllD(:,:,NoF+1:2*NoF)=AllDe;

Nvol=length(bval);
S=zeros(Nvol,1);
for s=1:Nvol
    Ss=0;
    ni=squeeze(bvec(:,s));
    for c=1:NoF*2
        adci=squeeze(AllD(:,:,c));
        adci=ni'*adci*ni;
        Ss=Ss+S0*(f(c)*exp(-bval(s)*adci));
    end
    S(s)=Ss;
end
if sig~=0
    %reset(RandStream.getGlobalStream,sum(100*clock))
    rng('shuffle')
    GR=sig*randn(Nvol,1);
    GI=sig*randn(Nvol,1);
    S=sqrt((S/sqrt(2)+GR).^2+(S/sqrt(2)+GI).^2);
end

%DT=AllD(:,:,1)*f(1)+AllD(:,:,2)*f(2)+AllD(:,:,3)*f(3)+AllD(:,:,4)*f(4);
%W=General_Lazar_Eq15_all_rnh(AllD,DT,f);
%[Dele,Wele,S0]=DKI_ULLS_fit(S,bval,bvec);

%A=[W',Wele]

