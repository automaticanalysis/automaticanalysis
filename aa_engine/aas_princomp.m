% Like princomp command in stats toolbox
% Similar output, except that it doesn't reorder components
%
function [coeff score latent]=aas_princomp(x)

% Zero columns, important to make SVD=princomp
x=x-repmat(mean(x,1),[size(x,1) 1]);

% SVD
[u s v]=svd(x);

% GET PCA
coeff=v;
score=u*s;
latent=diag(s).^2/(size(x,1)-1);
