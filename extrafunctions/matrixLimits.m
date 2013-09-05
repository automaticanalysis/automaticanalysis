% Find out an "inclusive mask" matrix's limits and save them to a cell array
function matrixLimit = matrixLimits(matrix, padding)
if nargin < 2
    padding = 0;
end

matrixSize = size(matrix);
dims =  length(matrixSize);

if length(padding) == 1
    padding = repmat(padding, [1 dims]);
end

matrixLimit = cell(dims,1);

for m = 1:dims
    matrixLimit{m} = matrix~=0 & isfinite(matrix);
    % For each dimension check which dimensions whose sum > 0
    for d = dims:-1:1
        if d ~= m
            matrixLimit{m} = squeeze(sum(matrixLimit{m},d));
        end
    end
    % Find first and last element that is > 0
    matrixLimit{m} = [find(matrixLimit{m}>0, 1, 'first') - padding(m), ...
        find(matrixLimit{m}>0, 1, 'last') + padding(m)];
    if matrixLimit{m}(1) < 1
        matrixLimit{m}(1) = 1;
    end
    if matrixLimit{m}(2) > matrixSize(m)
        matrixLimit{m}(2) = matrixSize(m);
    end
end