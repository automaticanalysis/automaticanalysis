function F=DTImodel_we(x, A)

XDiso = [3e-3; 3e-3; 3e-3; 0; 0; 0; x(7)];
F = (x(8).*exp(A*x(1:7)) + (1-x(8)).*exp(A*XDiso));