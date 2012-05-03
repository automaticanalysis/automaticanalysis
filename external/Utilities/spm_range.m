function y = spm_range(x,dim)
% Computes the difference between the min and max of a vector. If you need
% to use it on a matrix, then you need to specify which dimension to
% operate on.
if nargin < 2
    y = max(x) - min(x);
else
    y = max(x,[],dim) - min(x,[],dim);
end
