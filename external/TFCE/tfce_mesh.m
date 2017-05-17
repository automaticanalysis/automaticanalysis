function tfce = tfce_mesh(faces, t, deltaT)
% FORMAT tfce = tfce_mesh(faces, t, deltaT)
%
% Christian Gaser
% $Id: tfce_mesh.m 78 2015-09-10 10:23:01Z gaser $

E = 1.0;
H = 2.0;

tfce = zeros(size(t));
t_max = max(abs(t(:)));

n_steps = ceil(t_max/deltaT);

%-Compute the (reduced) adjacency matrix
%--------------------------------------------------------------------------
A       = spm_mesh_adjacency(faces);
A       = A + speye(size(A));

for j = 1:n_steps

  thresh = j*deltaT;

  % positive t-values
  T   = (t >= thresh);
  ind = find(T);
  
  C = find_connected_component(A, T);
  C = C(ind);
    
  for i = 1:max(C)
    M = find(C == i);
    k = length(M);
    tfce(ind(M)) = tfce(ind(M)) + power(k,E)*power(thresh,H)*T(ind(M));
  end
  
  
  % negative t-values
  T   = (-t >= thresh);
  ind = find(T);
  
  C = find_connected_component(A, T);
  C = C(ind);
    
  for i = 1:max(C)
    M = find(C == i);
    k = length(M);
    tfce(ind(M)) = tfce(ind(M)) - power(k,E)*power(thresh,H)*T(ind(M));
  end
  
end

% make tfce values independend from deltaT
tfce = tfce*deltaT;


%--------------------------------------------------------------------------
function C = find_connected_component(A, T);
% find connected components 
% FORMAT C = find_connected_component(A,T)
% A        - a [nxn[ (reduced) adjacency matrix
% T        - a [nx1] data vector (using NaNs or logicals), n = #vertices
%
% C        - a [nx1] vector of cluster indices
%
% modified version from spm_mesh_clusters.m 5065 2012-11-16 20:00:21Z guillaume
%


%-Input parameters
%--------------------------------------------------------------------------
if ~islogical(T)
  T   = ~isnan(T);
end
  
A1 = A;
A1(~T,:) = [];
A1(:,~T) = [];

%-And perform Dulmage-Mendelsohn decomposition to find connected components
%--------------------------------------------------------------------------
[p,q,r] = dmperm(A1);
N       = diff(r);
CC      = zeros(size(A1,1),1);
for i = 1:length(r)-1
  CC(p(r(i):r(i+1)-1)) = i;
end
C       = NaN(numel(T),1);
C(T)    = CC;

%-Sort connected component labels according to their size
%--------------------------------------------------------------------------
[N,ni]  = sort(N(:), 1, 'descend');
[ni,ni] = sort(ni);
C(T)    = ni(C(T));
