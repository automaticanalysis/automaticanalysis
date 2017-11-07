function [Da,ADe,RDe,tortu]=EQ12EQ13EQ14EQ15_Fieremans(DTe,DTa)

Da=trace(DTa);%Eq 12

[~,Ldte]=eig(DTe);
dLdte=diag(Ldte);
sLdte=sort(dLdte);
l1dte=sLdte(3);
l2dte=sLdte(2);
l3dte=sLdte(1);
ADe=l1dte; %Eq 13
RDe=(l2dte+l3dte)/2; %Eq 14

tortu=ADe/RDe; %Eq 15
