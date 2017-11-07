function ODF=DKIODF(Dv,Kv,x,alfa,neg,U)
% Compute the apparent value of the DKI-ODF along the direction x
% given voxel's Dv and Kv
% Inputs:
% Dv=[D11 D22 D33 D12 D13 D23];
% Kv=[W1111 W2222 W3333 W1112 W1113...
%     W1222 W2223 W1333 W2333 W1122...
%     W1133 W2233 W1123 W1223 W1233];
% Rafael N H
% 23/04/2014
% Update 06/04/2015 

if length(x)==2
    
    az=x(1);
    el=x(2);
    [xxx, yyy, zzz]=sph2cart(az,el,1);
    gnz=[xxx,yyy,zzz];
    
else
    
    gnz=[x(1),x(2),x(3)];
    
end

if nargin<5
    
    MD=mean(Dv(1:3)); % Mean diffusivity
    
    DT=[Dv(1),Dv(4),Dv(5);...
        Dv(4),Dv(2),Dv(6);...
        Dv(5),Dv(6),Dv(3) ];
    
    U=pinv(DT)*MD;% dimensionless tensor U
    
    neg=false;
    
elseif nargin == 5
    
    MD=mean(Dv(1:3)); % Mean diffusivity
    
    DT=[Dv(1),Dv(4),Dv(5);...
        Dv(4),Dv(2),Dv(6);...
        Dv(5),Dv(6),Dv(3) ];
    
    U=pinv(DT)*MD;% dimensionless tensor U
    
end

Un=U*gnz';
nUn=gnz*U*gnz';

GODF=(1/(nUn)).^((alfa+1)/2); % Gaussian ODF

% matrice to compute the DKIODF
V11= Un(1)^2/(nUn);
V22= Un(2)^2/(nUn);
V33= Un(3)^2/(nUn);
V12= Un(1)*Un(2)/(nUn);
V13= Un(1)*Un(3)/(nUn);
V23= Un(2)*Un(3)/(nUn);

ODF=GODF*(1+1/24*(Kv(1)*(3*U(1,1)*U(1,1)-6*(alfa+1)*U(1,1)*V11+(alfa+1)*(alfa+3)*V11*V11)+...
    Kv(2)*(3*U(2,2)*U(2,2)-6*(alfa+1)*U(2,2)*V22+(alfa+1)*(alfa+3)*V22*V22)+...
    Kv(3)*(3*U(3,3)*U(3,3)-6*(alfa+1)*U(3,3)*V33+(alfa+1)*(alfa+3)*V33*V33)+...
    Kv(4)*(12*U(1,1)*U(1,2)-12*(alfa+1)*U(1,1)*V12-12*(alfa+1)*U(1,2)*V11+4*(alfa+1)*(alfa+3)*V11*V12)+...
    Kv(5)*(12*U(1,1)*U(1,3)-12*(alfa+1)*U(1,1)*V13-12*(alfa+1)*U(1,3)*V11+4*(alfa+1)*(alfa+3)*V11*V13)+...
    Kv(6)*(12*U(1,2)*U(2,2)-12*(alfa+1)*U(1,2)*V22-12*(alfa+1)*U(2,2)*V12+4*(alfa+1)*(alfa+3)*V12*V22)+...
    Kv(7)*(12*U(2,2)*U(2,3)-12*(alfa+1)*U(2,2)*V23-12*(alfa+1)*U(2,3)*V22+4*(alfa+1)*(alfa+3)*V22*V23)+...
    Kv(8)*(12*U(1,3)*U(3,3)-12*(alfa+1)*U(1,3)*V33-12*(alfa+1)*U(3,3)*V13+4*(alfa+1)*(alfa+3)*V13*V33)+...
    Kv(9)*(12*U(2,3)*U(3,3)-12*(alfa+1)*U(2,3)*V33-12*(alfa+1)*U(3,3)*V23+4*(alfa+1)*(alfa+3)*V23*V33)+...
    Kv(10)*(6*U(1,1)*U(2,2)+12*U(1,2)*U(1,2)-6*(alfa+1)*U(1,1)*V22-6*(alfa+1)*U(2,2)*V11-24*(alfa+1)*U(1,2)*V12+2*(alfa+1)*(alfa+3)*V11*V22+4*(alfa+1)*(alfa+3)*V12*V12)+...
    Kv(11)*(6*U(1,1)*U(3,3)+12*U(1,3)*U(1,3)-6*(alfa+1)*U(1,1)*V33-6*(alfa+1)*U(3,3)*V11-24*(alfa+1)*U(1,3)*V13+2*(alfa+1)*(alfa+3)*V11*V33+4*(alfa+1)*(alfa+3)*V13*V13)+...
    Kv(12)*(6*U(2,2)*U(3,3)+12*U(2,3)*U(2,3)-6*(alfa+1)*U(2,2)*V33-6*(alfa+1)*U(3,3)*V22-24*(alfa+1)*U(2,3)*V23+2*(alfa+1)*(alfa+3)*V22*V33+4*(alfa+1)*(alfa+3)*V23*V23)+...
    Kv(13)*(12*U(1,1)*U(2,3)+24*U(1,2)*U(1,3)-12*(alfa+1)*U(1,1)*V23-12*(alfa+1)*U(2,3)*V11-24*(alfa+1)*U(1,2)*V13-24*(alfa+1)*U(1,3)*V12+4*(alfa+1)*(alfa+3)*V11*V23+8*(alfa+1)*(alfa+3)*V12*V13)+...
    Kv(14)*(12*U(2,2)*U(1,3)+24*U(2,1)*U(2,3)-12*(alfa+1)*U(2,2)*V13-12*(alfa+1)*U(1,3)*V22-24*(alfa+1)*U(1,2)*V23-24*(alfa+1)*U(2,3)*V12+4*(alfa+1)*(alfa+3)*V22*V13+8*(alfa+1)*(alfa+3)*V12*V23)+...
    Kv(15)*(12*U(3,3)*U(1,2)+24*U(3,1)*U(3,2)-12*(alfa+1)*U(3,3)*V12-12*(alfa+1)*U(1,2)*V33-24*(alfa+1)*U(1,3)*V23-24*(alfa+1)*U(2,3)*V13+4*(alfa+1)*(alfa+3)*V33*V12+8*(alfa+1)*(alfa+3)*V13*V23)));

if neg
    ODF=-ODF;
end