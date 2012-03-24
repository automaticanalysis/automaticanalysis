% Function for rescaling a file for coregistration...

function rescale4coreg(Cimg)

    V = spm_vol(Cimg);
    Y = spm_read_vols(V);
    
    % Anything that is less than 0 in original image gets set to 0...
    % Typically this is erroneous, but let's check how many voxels there are...
    if sum(Y < 0) > 0
        fprintf('\nThere are %d voxels below 0 in your dataset', sum(Y < 0))
        if sum(Y < 0) > 0.00001 * V.dim(1) * V.dim(2) * V.dim(3)
            error('That is an awful lot of voxels below 0')
        end
    end
    Y(Y <= 0) = 0;
    
    % Now, let's round the values from the image, since we don't care about the
    % decimal places in any case...
    % Check that the rounding is in a good direction...
    Y = round(Y);
    
    % Get the values that Y can take:
    Yvals = 1:1:max(Y(:));
    
    % Now, let's get the histogram of the values above 0...
    H = hist(Y(Y>0), Yvals);
        
    % Collapsed Y (where Y takes values 0:255)
    cY = Y;
    
    convertN = cell(1,255);
    
    maxN = 0;
    for n = 1:255
        booBin = 0;
        minN = maxN+1;
        maxN = minN; % At the start...
        while booBin == 0
            maxN = maxN + 1;
            
            binSum = sum(H(minN:end)) / (255 + 1 - n);
            % Once we get past the point... choose whether we do left or right of threshold...
            if sum(H(minN:maxN)) >= binSum
                if abs(sum(H(minN:maxN)) - binSum) > abs(sum(H(minN:(maxN-1))) - binSum)
                    maxN = maxN - 1;
                end
                convertN{n} = minN:maxN;
                cY(Y>=minN & Y<=maxN) = n;
                booBin = 1;
            end
        end
    end
    
    spm_write_vol(V,cY);
    
    fprintf('\tThe correlation between the original and rescaled image is... %0.2f\n', corr(Y(:), cY(:)));