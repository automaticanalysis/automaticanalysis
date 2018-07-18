function [MD,FA,l1,l2,l3,v1,v2,v3]=DTImetrics(DT)

[V,L]=eig(DT);
dL=diag(L);
[sL,is]=sort(dL);
l1=sL(3);
l2=sL(2);
l3=sL(1);

v1=V(:,is(3));
v2=V(:,is(2));
v3=V(:,is(1));

MD=(l1+l2+l3)/3;
FA=sqrt(3/2*((l1-MD)^2+(l2-MD)^2+(l3-MD)^2 )/(l1^2+l2^2+l3^2) );
