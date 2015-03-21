function [AK,RK]=fun_AK_RK(Dv,Kv,Direc)
% Computing AK and RK based on the mean direction of the DT (i.e.
% principal eigenvector)
% Implemented by Rafael Neto Henriques (Last review 20/02/2015)

N=size(Direc,1);
nodir=30; % parameter ajusted based on some visual inspections, corresponds 
% to the number of points to compute the integral of the prependicular 
% kurtosis
AK=zeros(N,1);
RK=AK;
MD=mean(Dv(1:3));
for vt=1:N
    
    % axial K
    gnz=Direc(vt,:);
    Dapp=gnz(1)*gnz(1)*Dv(1)+gnz(2)*gnz(2)*Dv(2)+gnz(3)*gnz(3)*Dv(3)+...
        2*gnz(1)*gnz(2)*Dv(4)+2*gnz(1)*gnz(3)*Dv(5)+2*gnz(2)*gnz(3)*Dv(6);
    
    Wapp=...
        gnz(1)*gnz(1)*gnz(1)*gnz(1)*Kv(1)+...
        gnz(2)*gnz(2)*gnz(2)*gnz(2)*Kv(2)+...
        gnz(3)*gnz(3)*gnz(3)*gnz(3)*Kv(3)+...
        4*gnz(1)*gnz(1)*gnz(1)*gnz(2)*Kv(4)+...
        4*gnz(1)*gnz(1)*gnz(1)*gnz(3)*Kv(5)+...
        4*gnz(1)*gnz(2)*gnz(2)*gnz(2)*Kv(6)+...
        4*gnz(2)*gnz(2)*gnz(2)*gnz(3)*Kv(7)+...
        4*gnz(1)*gnz(3)*gnz(3)*gnz(3)*Kv(8)+...
        4*gnz(2)*gnz(3)*gnz(3)*gnz(3)*Kv(9)+...
        6*gnz(1)*gnz(1)*gnz(2)*gnz(2)*Kv(10)+...
        6*gnz(1)*gnz(1)*gnz(3)*gnz(3)*Kv(11)+...
        6*gnz(2)*gnz(2)*gnz(3)*gnz(3)*Kv(12)+...
        12*gnz(1)*gnz(1)*gnz(2)*gnz(3)*Kv(13)+...
        12*gnz(1)*gnz(2)*gnz(2)*gnz(3)*Kv(14)+...
        12*gnz(1)*gnz(2)*gnz(3)*gnz(3)*Kv(15);
    
    Kapp=(MD/Dapp)^2*Wapp;
    AK(vt)=Kapp;
    
    % radial k
    if(Direc(vt,1)==1);
        Dir_per=circulo_def_direccao_y(Direc(vt,1),Direc(vt,2),Direc(vt,3),nodir);
    else
        Dir_per=circulo_def_direccao_x(Direc(vt,1),Direc(vt,2),Direc(vt,3),nodir);
    end
    
    rki=0;
    for s=1:nodir
        gnz=Dir_per(s,:);
        Dapp=gnz(1)*gnz(1)*Dv(1)+gnz(2)*gnz(2)*Dv(2)+gnz(3)*gnz(3)*Dv(3)+...
            2*gnz(1)*gnz(2)*Dv(4)+2*gnz(1)*gnz(3)*Dv(5)+2*gnz(2)*gnz(3)*Dv(6);
        
        Wapp=...
            gnz(1)*gnz(1)*gnz(1)*gnz(1)*Kv(1)+...
            gnz(2)*gnz(2)*gnz(2)*gnz(2)*Kv(2)+...
            gnz(3)*gnz(3)*gnz(3)*gnz(3)*Kv(3)+...
            4*gnz(1)*gnz(1)*gnz(1)*gnz(2)*Kv(4)+...
            4*gnz(1)*gnz(1)*gnz(1)*gnz(3)*Kv(5)+...
            4*gnz(1)*gnz(2)*gnz(2)*gnz(2)*Kv(6)+...
            4*gnz(2)*gnz(2)*gnz(2)*gnz(3)*Kv(7)+...
            4*gnz(1)*gnz(3)*gnz(3)*gnz(3)*Kv(8)+...
            4*gnz(2)*gnz(3)*gnz(3)*gnz(3)*Kv(9)+...
            6*gnz(1)*gnz(1)*gnz(2)*gnz(2)*Kv(10)+...
            6*gnz(1)*gnz(1)*gnz(3)*gnz(3)*Kv(11)+...
            6*gnz(2)*gnz(2)*gnz(3)*gnz(3)*Kv(12)+...
            12*gnz(1)*gnz(1)*gnz(2)*gnz(3)*Kv(13)+...
            12*gnz(1)*gnz(2)*gnz(2)*gnz(3)*Kv(14)+...
            12*gnz(1)*gnz(2)*gnz(3)*gnz(3)*Kv(15);
        
        Kapp=(MD/Dapp)^2*Wapp;
        rki=rki+Kapp;
    end
    RK(vt)=rki/nodir;
end
end

function Dir_samples=circulo_def_direccao_x(a,b,c,nodir)
%% Rafael Neto Henriques - calculo do circulo unitario, prependicular a uma
% direcao dada pelo utilizador (parametros de entreda a, b,c - x,y,z do 
% vector 
%%
% isto ainda pode ser simplificado por uma expressso unica, contudo para
% ser mais facil de detectar possiveis problemas resolvi deixar em passos

%[a,b,c]=sph2cart(theta,phi,1);
%ff=figure;
%plot3([0 a],[0 b],[0 c],'black')
alfa=pi/nodir:pi/nodir:pi;% half of points because is symetric
calfa=cos(alfa);
salfa=sin(alfa);
sq=sqrt(b^2+c^2);
R1=[0, -c/sq, b/sq];
R2=[-sq,b*a/sq,c*a/sq];
Dir_samples=zeros(nodir,3);
for s=1:nodir
    Dir_samples(s,:)=R1*(calfa(s))+R2*(salfa(s));
end
end

function Dir_samples=circulo_def_direccao_y(a,b,c,nodir)
%% Rafael Neto Henriques - calculo do circulo unitario, prependicular a uma
% direcao dada pelo utilizador (parametros de entreda a, b,c - x,y,z do 
% vector 
%%
% isto ainda pode ser simplificado por uma expressso unica, contudo para
% ser mais facil de detectar possiveis problemas resolvi deixar em passos

%[a,b,c]=sph2cart(theta,phi,1);
%ff=figure;
%plot3([0 a],[0 b],[0 c],'black')
alfa=pi/nodir:pi/nodir:pi;% half of points because is symetric
calfa=cos(alfa);
salfa=sin(alfa);
sq=sqrt(a^2+c^2);
R1=[-c/sq, 0, a/sq];
R2=[-a*b/sq,sq,-c*b/sq];
Dir_samples=zeros(nodir,3);
for s=1:nodir
    Dir_samples(s,:)=R1*(calfa(s))+R2*(salfa(s));
end
end