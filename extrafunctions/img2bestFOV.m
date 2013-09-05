function coords = img2bestFOV(Mimg, sizeFOV, binarise)
if nargin == 2
    binarise = 0;
elseif nargin < 2
    error('Not enough input agruments')
end

% Load image first
V = spm_vol(Mimg);
M = spm_read_vols(V);

if binarise
    % Round the image...
    M = round(M);
end

for dim = 1:3
    tM = M;
    
    for d = 1:3
        if d ~= dim
            tM = mean(tM,d);
        end
    end
    tM = squeeze(tM);
    
    % DEBUG
    %{
    try close(3); catch; end
    figure(3)
    plot(tM)
    %}
    
    sumM = zeros(1, (V.dim(dim) - sizeFOV(dim)));
    if ~isempty(sumM)
    for s = 1:length(sumM)
        sumM(s) = sum(tM(s:(s+sizeFOV(dim)-1)));
    end   
        coords(dim) = find(max(sumM)==sumM, 1) + sizeFOV(dim) / 2;
    else
        coords(dim) = V.dim(dim) ./ 2;
    end
end