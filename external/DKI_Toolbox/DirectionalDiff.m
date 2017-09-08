function Dapp=DirectionalDiff(Dv,x,neg)
% Compute the value of diffusion along the direction x
% given voxel's Diffusion tensor Dv
% Dv=[D11 D22 D33 D12 D13 D23];

% Rafael N H
% 07/04/2015

if nargin < 3
    neg=false;
end

if length(x)==2
    az=x(1);
    el=x(2);
    [xxx, yyy, zzz]=sph2cart(az,el,1);
    gnz=[xxx,yyy,zzz];
else 
    gnz=[x(1),x(2),x(3)]; 
end

Dapp=gnz(1)*gnz(1)*Dv(1)+gnz(2)*gnz(2)*Dv(2)+gnz(3)*gnz(3)*Dv(3)+...
    2*gnz(1)*gnz(2)*Dv(4)+2*gnz(1)*gnz(3)*Dv(5)+2*gnz(2)*gnz(3)*Dv(6);

if neg
    Dapp=-Dapp;
end