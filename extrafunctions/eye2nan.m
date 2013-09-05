% Supplants leading diagonal of square matrix with NaNs
function M = eye2nan(M)
sM = size(M);
if length(sM) ~= 2 || sM(1) ~= sM(2)
   error('Not a square matrix') 
end

E = logical(eye(size(M)));
M(E) = NaN;