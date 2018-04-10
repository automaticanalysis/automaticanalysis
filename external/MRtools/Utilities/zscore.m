function [z,mu,sigma] = zscore(x,flag,dim)
if isequal(x,[]), z = []; return; end

if nargin < 2
    flag = 0;
end
if nargin < 3
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

mu = mean(x,dim);
sigma = std(x,flag,dim);
sigma0 = sigma;
sigma0(sigma0==0) = 1;
z = bsxfun(@minus,x, mu);
z = bsxfun(@rdivide, z, sigma0);

