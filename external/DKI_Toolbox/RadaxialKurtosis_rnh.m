function [AK,RK]=RadaxialKurtosis_rnh(Dv,Kv,Direc)

N=size(Direc,1);
AK=zeros(N,1);
RK=AK;
MD=mean(Dv(1:3));
for vt=1:N
    
    % axial K
    gnz=Direc(vt,:);
    Kapp=DirectionalKurt(Dv,Kv,gnz,false,MD);
    AK(vt)=Kapp;
    
    % radial K
    Kapp=RadialKurtosis(Dv,Kv,gnz,false,MD);
    RK(vt)=Kapp;
end