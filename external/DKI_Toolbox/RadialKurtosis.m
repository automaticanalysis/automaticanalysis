function RK=RadialKurtosis(Dv,Kv,Direc,neg,MD,nodir,calfa,salfa)

if nargin <6
    nodir=30; 
    alfa=pi/nodir:pi/nodir:pi;% half of points because is symetric
    calfa=cos(alfa)';
    salfa=sin(alfa)';
    if nargin <5
        MD=mean(Dv(1:3));
        if naring <4
            neg=false;
        end
    end
end

if size(Direc,2)==2
    az=Direc(1);
    el=Direc(2);
    [xxx, yyy, zzz]=sph2cart(az,el,1);
    gnz=[xxx,yyy,zzz];
else 
    gnz=[Direc(1),Direc(2),Direc(3)]; 
end


% radial k
if(Direc(1,1)==1);
    Dir_per=Perp_cirle_y(gnz(1),gnz(2),gnz(3),nodir,calfa,salfa);
else
    Dir_per=Perp_circle_x(gnz(1),gnz(2),gnz(3),nodir,calfa,salfa);
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

RK=rki/nodir;

if neg
    RK=-RK;
end
