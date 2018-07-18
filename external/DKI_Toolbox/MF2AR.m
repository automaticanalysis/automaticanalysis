function [A, R]=MF2AR(M,F)

A=M*(1+(2*F)/sqrt(3-2*F^2));
R=(3*M-A)/2;

end