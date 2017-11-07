function ODF=DKIODF_v(Dv,Kv,V,alfa)
% Compute the apparent value of the DKI-ODF along the directions 
% in matrix V given voxel's Dv and Kv
% Inputs:
% Dv=[D11 D22 D33 D12 D13 D23];
% Kv=[W1111 W2222 W3333 W1112 W1113...
%     W1222 W2223 W1333 W2333 W1122...
%     W1133 W2233 W1123 W1223 W1233];
% V matrix Nx3 with the coordinates x, y, z of N directions
% Rafael N H
% 23/04/2014
% Update 06/04/2015 (Unifying functions)

N=size(V,1);

ODF=zeros(N,1); % values of DKI ODF

MD=mean(Dv(1:3)); % Mean diffusivity

DT=[Dv(1),Dv(4),Dv(5);...
    Dv(4),Dv(2),Dv(6);...
    Dv(5),Dv(6),Dv(3) ];

U=pinv(DT)*MD;% dimensionless tensor U

for vt=1:N
    ODF(vt)=DKIODF(Dv,Kv,V(vt,:),alfa,false,U);
end