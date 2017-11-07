function RK=DirectionalKurt_V(Dv,Kv,V)
% Compute the value of kurtosis along the directions 
% in matrix V given voxel's Dv and Kv
% Inputs:
% Dv=[D11 D22 D33 D12 D13 D23];
% Kv=[W1111 W2222 W3333 W1112 W1113...
%     W1222 W2223 W1333 W2333 W1122...
%     W1133 W2233 W1123 W1223 W1233];
% V matrix Nx3 with the coordinates x, y, z of N directions
% Rafael N H
% 07/04/2015

N=size(V,1);
RK=zeros(N,1);

for vt=1:N
    x=V(vt,:);
    Kapp=DirectionalKurt(Dv,Kv,x);
    RK(vt)=Kapp;
end



