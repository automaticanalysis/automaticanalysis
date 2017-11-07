function [MK, AK, RK]=DKImetrics(D,W,V)
% Computes standard DKI metrics from single voxel DT and KT
% Rafael Neto Henriques 04/04/2015

% Mean Kurtosis
MK=MeanKurtosis_rnh(D,W,V);

% radial e axial kurtosis
Dv=[D(1) D(4) D(5);...
    D(4) D(2) D(6);...
    D(5) D(6) D(3)];
[Vecs,L]=eig(Dv);
dL=diag(L);
[da,is]=sort(dL);
v1=Vecs(:,is(3));
[AK,RK]=RadaxialKurtosis_rnh(D,W,v1');
