function [MD, FA, AD, RD, l1, l2, l3, v1, v2, v3] = fun_DTI_metrics(DT,Mask)
% Compute standard DTI metrics from the DT
% Rafael Neto Henriques (last review 23/02/2015)

%Initialization of variables
[N1,N2,N3]=size(Mask);
ZER1=zeros(N1,N2,N3);
ZER2=zeros(N1,N2,N3,3);
MD=ZER1;
FA=ZER1;
RD=ZER1;
AD=ZER1;
l1=ZER1;
l2=ZER1;
l3=ZER1;
v1=ZER2;
v2=ZER2;
v3=ZER2;

for k=1:N3
    for j=1:N2
        for i=1:N1
            if(Mask(i,j,k)==1)
                D=DT(i,j,k,:);
                Dv=[D(1) D(4) D(5);...
                    D(4) D(2) D(6);...
                    D(5) D(6) D(3)];
                [MDi,FAi,l1i,l2i,l3i,v1i,v2i,v3i]=IVdiff(Dv);
                
                MD(i,j,k)=MDi;
                FA(i,j,k)=FAi;
                l1(i,j,k)=l1i;
                l2(i,j,k)=l2i;
                l3(i,j,k)=l3i;
                AD(i,j,k)=l1i;
                RD(i,j,k)=(l2i+l3i)/2;
                v1(i,j,k,:)=v1i;
                v2(i,j,k,:)=v2i;
                v3(i,j,k,:)=v3i;
            end
        end
    end
    %disp(k)
end
end


function [MD,FA,l1,l2,l3,v1,v2,v3]=IVdiff(DT)

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
end