function [M, F]=AR2MF(A,R)

M=(A+2*R)/3;

F=sqrt(3/2*((A-M)^2+2*(R-M)^2)/(A^2+2*R^2));

end