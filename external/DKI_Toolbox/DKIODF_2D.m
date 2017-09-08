function RODF=DKIODF_2D(Dv,Kv,theta,phi,alfa)
% Compute the apparent value of the DKI-ODF along the directions which coordinates are
% saved on the meshes theta and phi
%
% Inputs:
% Dv=[D11 D22 D33 D12 D13 D23];
% Kv=[W1111 W2222 W3333 W1112 W1113...
%     W1222 W2223 W1333 W2333 W1122...
%     W1133 W2233 W1123 W1223 W1233];
% Rafael N H
% 23/04/2014
% Update 07/04/2015 (Unifying functions)

[x,y,z]=sph2cart(theta,phi,1);
N=size(x);
na=N(1);
nb=N(2);

RODF=zeros(N); % values of DKI ODF

MD=mean(Dv(1:3)); % Mean diffusivity

DT=[Dv(1),Dv(4),Dv(5);...
    Dv(4),Dv(2),Dv(6);...
    Dv(5),Dv(6),Dv(3) ];

U=pinv(DT)*MD;% dimensionless tensor U

for a=1:na
    for b=1:nb
        gnz=[x(a,b),y(a,b),z(a,b)];   
        RODF(a,b)=DKIODF(Dv,Kv,gnz,alfa,false,U);
    end
end