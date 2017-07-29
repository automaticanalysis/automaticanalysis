function RK=DirectionalKurt_2D(Dv,Kv,theta,phi)
% Compute the value of diffusion along the directions which coordinates are
% saved on the meshes theta and phi
% 
% Inputs:
% Dv=[D11 D22 D33 D12 D13 D23];
% Kv=[W1111 W2222 W3333 W1112 W1113...
%     W1222 W2223 W1333 W2333 W1122...
%     W1133 W2233 W1123 W1223 W1233];
% 07/04/2015

[x,y,z]=sph2cart(theta,phi,1);
N=size(x);
na=N(1);
nb=N(2);

RK=zeros(N);

for a=1:na
    for b=1:nb
        gnz=[x(a,b),y(a,b),z(a,b)];        
        RK(a,b)=DirectionalKurt(Dv,Kv,gnz);
    end
end



