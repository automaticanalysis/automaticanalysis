function DKIu_simulations

% Ask number of fibres
ppp={'Number of fibers:'};
def={'2'};
r=inputdlg(ppp,'Imput required',1,def);
NoF=str2double(r{1}); %Number of Fibers

% Ask the values of intra and extra celular diffusivity
adbie={'AD (intra comp):','RD (intra comp):',...
    'AD (extra comp):','RD (extra comp):',...
    'f (intra vs extra):'};
def={'0.99e-3','0',...
    '2.26e-3','0.87e-3',...
    '0.49'};
r=inputdlg(adbie,'Intra and Extra celular parameters',1,def);
DAi=str2double(r{1});
DRi=str2double(r{2});
DAe=str2double(r{3});
DRe=str2double(r{4});
fie=str2double(r{5});
Di=[DAi 0 0;0 DRi 0;0 0 DRi];
De=[DAe 0 0;0 DRe 0;0 0 DRe];

% Inicilize parameter for plots
np=60;
theta=linspace(0,2*pi,np);
phi=linspace(-pi/2,pi/2,np);
[theta,phi]=meshgrid(theta,phi);

% Ask missing parameters for each fiber
def={num2str(1/NoF),'0','0'};

% inicialize variables
fp=zeros(1,NoF);
fphi=fp;
ftheta=fp;
AllDi=zeros(3,3,NoF);
AllDe=AllDi;
DTci=AllDi;
DTce=AllDi;
f_p1=fp; % compartment fractions part 1 (intracelular compartments)
f_p2=fp; % compartment fractions part 2 (extracelular compartments)

% parameters for fiber direction estimation procedures
load Dir125.mat
alfa=4;
MDc=0.8;
optionsT = optimset('TolX',1e-4); 

% Ask
for c=1:NoF
    if c==NoF
        rfp=1-sum(fp);
        if rfp<0
            errordlg('sum of fp for all compartments cannot be larger than 1')
            error('sum of f for all compartments cannot be larger than 1')
        else
            apf{1}='fiber f (use the recomended value):';
            def{1}=num2str(rfp);
        end
    end
    apf={['direction of fiber population #',num2str(c)],'phi (rad):','theta (rad):'};
    r=inputdlg(apf,['Fiber ', num2str(c)],1,def);
    
    % Fractions
    fp(1,c)=str2double(r{1});
    f_p1(1,c)=fp(1,c)*fie; 
    f_p2(1,c)=fp(1,c)*(1-fie);
    
    azi=str2num(r{2});
    ele=str2num(r{3});
    
    fphi(1,c)=azi;
    ftheta(1,c)=ele;
    
    [xa,ya,za]=sph2cart(azi,ele,1);
    Xa=[xa,ya,za];
    
    if azi == pi/2 && ele ==0
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
    
    Evector=[Xa',Xb',Xc'];
    
    Dci=Evector*Di*pinv(Evector);
    AllDi(:,:,c)=Dci;
    DTci(:,:,c)=f_p1(1,c)*Dci; % Matrix to later compute DT;
    
    Dce=Evector*De*pinv(Evector);
    AllDe(:,:,c)=Dce;
    DTce(:,:,c)=f_p2(1,c)*Dce; % Matrix to later compute DT;
end

% Compute DT, f and AllD
f=[f_p1 f_p2];
DTc=zeros(3,3,NoF*2);
AllD=DTc;
DTc(:,:,1:NoF)=DTci;
DTc(:,:,NoF+1:2*NoF)=DTce;
DT=sum(DTc,3);
AllD(:,:,1:NoF)=AllDi;
AllD(:,:,NoF+1:2*NoF)=AllDe;

% DT in vector form
Dv=[DT(1,1) DT(2,2) DT(3,3) DT(1,2) DT(1,3) DT(2,3)];

% Compute elements of KT
Kv=sim_kt(AllD,DT,f);

% plot do DT
r=DirectionalDiff_2D(Dv,theta,phi);
[xx,yy,zz]=sph2cart(theta,phi,r);
fig=figure;
set(fig,'color',[1 1 1])
surf(xx,yy,zz,r)
title('DT'), xlabel('x'),ylabel('y'),zlabel('z')
view(0,0)

[vvv,lll]=eig(DT);
[ada,idi]=max(diag(lll));
FDir1=vvv(:,idi);
hold on
DTmax=max(r(:));
DTmax=5/3*DTmax;
plot3([DTmax*-1*FDir1(1) DTmax*1*FDir1(1)],...
    [DTmax*-1*FDir1(2) DTmax*1*FDir1(2)],...
    [DTmax*-1*FDir1(3) DTmax*1*FDir1(3)],...
    'color',[0 0 1],'LineWidth',2)
for fd=1:NoF
    [fx,fy,fz]=sph2cart(fphi(fd),ftheta(fd),1);
    plot3([DTmax*-0.9*fx, DTmax*0.9*fx],...
        [DTmax*-0.9*fy, DTmax*0.9*fy],...
        [DTmax*-0.9*fz, DTmax*0.9*fz],...
        'color',[0 0 0],'LineWidth',2)
end
xlim([-DTmax DTmax])
ylim([-DTmax DTmax])
zlim([-DTmax DTmax])

% plot do WT
r=DirectionalKurt_2D(Dv,Kv,theta,phi);
[xx,yy,zz]=sph2cart(theta,phi,r);
fig=figure;
set(fig,'color',[1 1 1])
surf(xx,yy,zz,r)
title('KT'), xlabel('x'),ylabel('y'),zlabel('z')
view(0,0)

kmax=max(r(:));
kmax=5/3*kmax;
[MV,MR]=Procedure_KT(Dv,Kv,V,Nvis,Uvis,optionsT,MDc);
hold on
for fd=1:length(MR)
    plot3([kmax*-1*MV(fd,1) kmax*1*MV(fd,1)]*2/3,...
        [kmax*-1*MV(fd,2) kmax*1*MV(fd,2)]*2/3,...
        [kmax*-1*MV(fd,3) kmax*1*MV(fd,3)]*2/3,...
        'color',[0 0 1],'LineWidth',2)
end
for fd=1:NoF
    [fx,fy,fz]=sph2cart(fphi(fd),ftheta(fd),1);
    plot3([kmax*-0.9*fx, kmax*0.9*fx]*2/3,...
        [kmax*-0.9*fy, kmax*0.9*fy]*2/3,...
        [kmax*-0.9*fz, kmax*0.9*fz]*2/3,...
        'color',[0 0 0],'LineWidth',2)
end
xlim([-kmax kmax])
ylim([-kmax kmax])
zlim([-kmax kmax])

% DKI-ODF
r=DKIODF_2D(Dv,Kv,theta,phi,alfa);
[xx,yy,zz]=sph2cart(theta,phi,r);
fig=figure;
set(fig,'color',[1 1 1])
surf(xx,yy,zz,r)
title('DKI-ODF'), xlabel('x'),ylabel('y'),zlabel('z')
view(0,0)

odfmax=max(r(:));
odfmax=5/3*odfmax;
[MV,MR]=Procedure_DKIODF(Dv,Kv,V,Nvis,Uvis,alfa,optionsT);         
hold on
for fd=1:length(MR)
    plot3([odfmax*-1*MV(fd,1) odfmax*1*MV(fd,1)],...
        [odfmax*-1*MV(fd,2) odfmax*1*MV(fd,2)],...
        [odfmax*-1*MV(fd,3) odfmax*1*MV(fd,3)],...
        'color',[0 0 1],'LineWidth',2)
end
for fd=1:NoF
    [fx,fy,fz]=sph2cart(fphi(fd),ftheta(fd),1);
    plot3([odfmax*-0.9*fx, odfmax*0.9*fx],...
        [odfmax*-0.9*fy, odfmax*0.9*fy],...
        [odfmax*-0.9*fz, odfmax*0.9*fz],...
        'color',[0 0 0],'LineWidth',2)
end
xlim([-odfmax odfmax])
ylim([-odfmax odfmax])
zlim([-odfmax odfmax])

