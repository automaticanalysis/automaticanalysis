function [y]=aas_normcdf(x)
y=0.5*erfc(-x/sqrt(2));