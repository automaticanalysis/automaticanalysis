function Y=nansum(X,dim)
if nargin<2 || ~isnumeric(dim)
    dim=1;
end
if dim>numel(size(X))
    error('Summing dimension is larger than the dimensions of X')
end
X(isnan(X))=0;
Y=sum(X,dim);
end