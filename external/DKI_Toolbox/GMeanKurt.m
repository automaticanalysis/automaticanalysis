function Mk=GMeanKurt(DMAT,WMAT,b,g,Md,mode)
% caluclate mean kurtosis given the Difusion and Kurtosis tensor
% g is the gradient direction and Md the mean difusion

% mode 2 if is given the values of W*Md*Md in WMAT matrix instead of W only,
% useful to avoid divisions from zero
% Rafael Henriques


ind=(b~=0);
gnz=g(:,ind);
    
Nn=sum(ind);
    
Mk=0;
for v=1:Nn
    %MK
    Dapp=0;
    for i=1:3
        for j=1:3
            Dapp=Dapp+gnz(i,v)*gnz(j,v)*DMAT(i,j);
        end
    end
    
    Wapp=0;
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                    Wapp=Wapp+gnz(i,v)*gnz(j,v)*gnz(k,v)*gnz(l,v)*WMAT(i,j,k,l);
                end
            end
        end
    end
    
    Kapp=1/(Dapp*Dapp)*Wapp;
    Mk=Mk+Kapp;
       
end

Mk=Mk/Nn;

if mode~=2
    Mk=Mk*Md*Md;
end

    