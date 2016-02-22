function [MK, AK, RK, AWF, ADa, ADe, RDe, tortu] = fun_DKI_metrics(DT,KT,Mask)
% Rafael Neto Henriques, 10/04/2014 MRC-CBU

% Dir125
V = Dir125('V');

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

optionsT = optimset('TolX',1e-2,'Display', 'off');

for k=1:N3
    for j=1:N2
        for i=1:N1
            if(Mask(i,j,k)==1)
                D=DT(i,j,k,:);
                Dv=[D(1) D(4) D(5);...
                    D(4) D(2) D(6);...
                    D(5) D(6) D(3)];
                W=KT(i,j,k,:);
                [Vecs,L]=eig(Dv);
                dL=diag(L);
                [da,is]=sort(dL);
                v1=Vecs(:,is(3));
                % Mean Kurtosis
                MK(i,j,k)=fun_MK(D,W,V);
                % radial e axial kurtosis
                [AK(i,j,k),RK(i,j,k)]=fun_AK_RK(D,W,v1');
                
                % advanced metrics
                [AWF(i,j,k),Dai,Dei,tortu(i,j,k),ADa(i,j,k),ADe(i,j,k),RDe(i,j,k),Kmaxi]=run_DKI_single_fiber_model(D,W,V,Dir125('Nvis'),Dir125('Uvis'),optionsT);
            end
        end
    end
end

% Thresholds added in 27/03/2014 (to remove huge negative implausible
% outliers
MK(MK>12)=12;
MK(MK<-10)=-10;
AK(AK>12)=12;
AK(AK<-10)=-10;
RK(RK>12)=12;
RK(RK<-10)=-10;