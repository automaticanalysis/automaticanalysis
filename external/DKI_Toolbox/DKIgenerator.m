function [KT,DT]=DKIgenerator(sel_plot,DAi,DRi,DAe,DRe,fie,fp,fphi,ftheta)

% archstr = computer('arch');
% st_comp=archstr(1:3);
% if strcmp(st_comp,'win')
%     sc='\';
% else
%     sc='/';
% end
% 
% addpath('..',sc,'ToolDKI_Mac')

% Inicilize parameter for plots
if nargin==0
    sel_plot=false;
end

if sel_plot==true
    np=60;
    theta=linspace(0,2*pi,np);
    phi=linspace(-pi/2,pi/2,np);
    [theta,phi]=meshgrid(theta,phi);
end

NoF = length(fp);
%simulation
%% Ask number of compartments
if nargin<2
ppp={'Number of fibers:'};
def={'2'};
r=inputdlg(ppp,'Imput required',1,def);
NoF=str2double(r{1});

%Number of Fibers
% Ask diferences between the intra and extra celular compartments
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
end

Di=[DAi 0 0;0 DRi 0;0 0 DRi];
De=[DAe 0 0;0 DRe 0;0 0 DRe];

if nargin<2
% Ask missing parameters for each fiber
apf={'fiber f','fiber phi (rad):','fiber theta (rad):'};
def={num2str(1/NoF),'0','0'};
fphi=fp;
ftheta=fp;
end

% inicialize variables
% fp=zeros(1,NoF);
AllDi=zeros(3,3,NoF);
AllDe=AllDi;
DTci=AllDi;
DTce=AllDi;
f_p1=fp; % compartment fractions part 1 (intracelular compartments)
f_p2=fp; % compartment fractions part 2 (extracelular compartments)


for c=1:NoF
    if nargin<2
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
    r=inputdlg(apf,['Fiber ', num2str(c)],1,def);
    
    % Fractions
    fp(1,c)=str2double(r{1});
    f_p1(1,c)=fp(1,c)*fie; % Hurray probabilities!!!
    f_p2(1,c)=fp(1,c)*(1-fie);    
    azi=str2num(r{2});
    ele=str2num(r{3});    
    fphi(1,c)=azi;
    ftheta(1,c)=ele;
    else 
        f_p1(c)=fp(c)*fie;
        f_p2(c)=fp(c)*(1-fie);
    end

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

% matrix para plot posterior
Dv=[DT(1,1) DT(2,2) DT(3,3) DT(1,2) DT(1,3) DT(2,3)];

% Calculo do tensor de kurtosis
W=General_Lazar_Eq15_all_rnh(AllD,DT,f);
%WT=Wcons(W);

if sel_plot==true
% plot do DT
r=IMPsdif3D(Dv,theta,phi);
[xx,yy,zz]=sph2cart(theta,phi,r);
fig=figure;
set(fig,'color',[1 1 1])
surf(xx,yy,zz,r)
title('DT'), xlabel('x'),ylabel('y'),zlabel('z')
view(0,0)

% plot do WT
r=IMPskur3D(Dv,W,theta,phi);
[xx,yy,zz]=sph2cart(theta,phi,r);
fig=figure;
set(fig,'color',[1 1 1])
surf(xx,yy,zz,r)
title('WT'), xlabel('x'),ylabel('y'),zlabel('z')
view(0,0)
end 

%fit

%yn=observations;

%Z=S(bn,gn,X)*yn;
%J = besselj(0,Z); 
DT=Dv;
KT=W;