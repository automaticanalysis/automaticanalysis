function RK=RadialKurtosis_v(Dv,Kv,Direc,MD,nodir,calfa,salfa)

if nargin <5
    nodir=30; 
    alfa=pi/nodir:pi/nodir:pi;% half of points because is symetric
    calfa=cos(alfa)';
    salfa=sin(alfa)';
    if nargin <4
        MD=mean(Dv(1:3));
    end
end

N=size(Direc,1);
RK=zeros(N,1);
for vt=1:N
    % radial k
    Kapp=RadialKurtosis(Dv,Kv,Direc(vt,:),false,MD,nodir,calfa,salfa);
    RK(vt)=Kapp;
end