function Wijkl_i=Wijkl(AllD,DT,f,i,j,k,l)
% Note Wijkl have to be divided by MD^2 to be values of kurtosis
% This function is the fraction with the same weight

ncomp=size(AllD,3);
Wijkl_fraq=0;
for comp=1:ncomp
    Dc=squeeze(AllD(:,:,comp));
    Wijkl_fraq=Wijkl_fraq+f(1,comp)*(Dc(i,j)*Dc(k,l)+Dc(i,k)*Dc(j,l)+Dc(i,l)*Dc(j,k));
end
Wijkl_i=Wijkl_fraq-DT(i,j)*DT(k,l)-DT(i,k)*DT(j,l)-DT(i,l)*DT(j,k);
