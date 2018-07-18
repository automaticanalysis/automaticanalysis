function [MK, AK, RK, AWF, ADa, ADe, RDe, tortu] = fun_DKI_metrics(DT,KT,Mask)
% Rafael Neto Henriques, 10/04/2014 MRC-CBU

%Initialization of variables
[N1,N2,N3]=size(Mask);
ZER1=zeros(N1,N2,N3);
MK=ZER1;
RK=ZER1;
AK=ZER1;
AWF=ZER1;
ADa=ZER1;
ADe=ZER1;
RDe=ZER1;
tortu=ZER1;


load('Dir125.mat')
optionsT = optimset('TolX',1e-2,'Display', 'off');

for k=1:N3
    for j=1:N2
        for i=1:N1
            if(Mask(i,j,k)==1)
                D=DT(i,j,k,:);
                W=KT(i,j,k,:);
                [MKi, AKi, RKi]=DKImetrics(D,W,V);
                MK(i,j,k)=MKi;
                RK(i,j,k)=RKi;
                AK(i,j,k)=AKi;
                
                % advanced metrics
                [AWF(i,j,k),Dai,Dei,tortu(i,j,k),ADa(i,j,k),ADe(i,j,k),RDe(i,j,k),Kmaxi]=DKI_single_fiber_model_rnh(D,W,V,Nvis,Uvis,optionsT);
            end
        end
    end
    disp(k)
end

% Thresholds added in 27/03/2014
MK(MK>12)=12;
MK(MK<-10)=-10;
AK(AK>12)=12;
AK(AK<-10)=-10;
RK(RK>12)=12;
RK(RK<-10)=-10;