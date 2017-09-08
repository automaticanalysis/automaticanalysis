function [MD, FA, AD, RD, l1, l2, l3, v1, v2, v3, CL, CP, CS, SVP] = fun_DTI_metrics(DT,Mask)

%inicializacao de variaveis
[N1,N2,N3, daaa]=size(Mask);
ZER1=zeros(N1,N2,N3);
ZER2=zeros(N1,N2,N3,3);
MD=ZER1;
FA=ZER1;
l1=ZER1;
l2=ZER1;
l3=ZER1;
RD=ZER1;
v1=ZER2;
v2=ZER2;
v3=ZER2;

for k=1:N3
    for j=1:N2
        for i=1:N1
            if(Mask(i,j,k)==1)
                D = DT(i,j,k,:);
                Dv=[D(1) D(4) D(5);...
                    D(4) D(2) D(6);...
                    D(5) D(6) D(3)];
                [MDi,FAi,l1i,l2i,l3i,v1i,v2i,v3i]=DTImetrics(Dv);
                
                MD(i,j,k)=MDi;
                FA(i,j,k)=FAi;
                l1(i,j,k)=l1i;
                l2(i,j,k)=l2i;
                l3(i,j,k)=l3i;
                RD(i,j,k)=(l2i+l3i)/2;
                v1(i,j,k,:)=v1i;
                v2(i,j,k,:)=v2i;
                v3(i,j,k,:)=v3i;
            end
        end
    end
    disp(k)
end

AD = l1;

%% Extra invariant maps
CL=(AD-l2)./AD;
CP=(l2-l3)./AD;
CS=l3./AD;

CP(isnan(CP))=0;
CL(isnan(CL))=0;
CS(isnan(CS))=0;

SVP=((CL>=0.4)&(CP<=0.2)&(CS<=0.35));
