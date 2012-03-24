function [p t] = corr2pt(C, N)

t = C ./ sqrt((1-C.^2) ./ (N-2));
p = 1 - tcdf(abs(t), N-2);