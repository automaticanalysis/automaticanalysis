function ROIdata = roi2hist(dataImg, ROIimg, threshROI, legendExtra)
if nargin < 2
   aas_log(aap,true, 'Requires at least 2 outputs...') 
end
if ischar(ROIimg)
    ROIimg = strvcat2cell(ROIimg);
end
if nargin < 3
    threshROI = [];
end
if nargin < 4 || isempty(legendExtra)
   legendExtra =  cell(1,numel(ROIimg));
end

% Load the data image
if isstr(dataImg)
    Y = spm_read_vols(spm_vol(dataImg));
else
    Y = dataImg;
end

ROIname = cell(1,numel(ROIimg));
ROIvol = cell(1,numel(ROIimg));
ROIdata = cell(1,numel(ROIimg));

for r = 1:length(ROIimg)    
    [junk, ROIname{r}] = fileparts(ROIimg{r});
    
    % Now load each of the ROIs we wish to examine (usually using the grey matter)
    rV = spm_vol(ROIimg{r});
    ROIvol{r} = spm_read_vols(rV);
    if isempty(threshROI)
        ROIvol{r} = round(ROIvol{r});
    else
        ROIvol{r} = ROIvol{r} > threshROI;
    end
    if length(unique(ROIvol{r})) > 2
        aas_log([],true,'Your ROI is not binary, set a threshold!')
    end
    
    if any(size(ROIvol{r})~=size(Y))
        aas_log([],true,'The dimensions of the data and the ROI do not match')
    end
end

for r = 1:numel(ROIimg)
    % Now get the voxels specific to each ROI
    ROIdata{r} = Y(ROIvol{r});
    ROIdata{r} = ROIdata{r}(isfinite(ROIdata{r}) & ROIdata{r} ~= 0); % We don't want to include zero values...
end

%% tSNR results figure!
colorsB = aas_colours;

% We need to make a string for legend...
legStr = cell(1,numel(ROIimg));
for r = 1:numel(ROIimg)
    legStr{r} = sprintf('%s; mn=%.2f; SD=%.2f; med=%.0f; (%.0fv) %s', ...
        ROIname{r},  ...
        mean(ROIdata{r}), ...
        std(ROIdata{r}), ...
        median(ROIdata{r}), ...
        length(ROIdata{r}), ...
        legendExtra{r});
end

try close(2); catch; end

figure(2)
set(2, 'Position', [0 0 1000 550])
minI = Inf;
maxI = -Inf;
windI = 0;
maxV = 0;
hold on

for r = 1:numel(ROIimg)
    % What range do the SNR values take?
    maxI = max(max(ROIdata{r}), maxI);
    minI = min(min(ROIdata{r}), minI);
    % What window do we wish to present?
    windI = max(median(ROIdata{r}) + std(ROIdata{r}) * 3, windI);
end
vals = floor(minI):windI/250:ceil(maxI);
for r = 1:numel(ROIimg)
    % Now make a histogram and "normalise" it
    H = hist(ROIdata{r}, vals);
    H = H./sum(H);
    % And decide what is the greatest prop value
    maxV = max(max(H), maxV);
    
    % Make bars semi-transparent for overlaying
    B = bar(vals, H, 1, 'FaceColor', colorsB{r});
    ch = get(B,'child');
    set(ch, 'faceA', 0.3, 'edgeA', 0.2);
end
vals = floor(minI):windI/100:ceil(windI);
% Set the axis to a good value!
axis([vals(1), vals(end), 0, maxV*1.05])
xlabel('Voxel value')
ylabel('Proportion of voxels')
set(gca,'XTick', 0:ceil(maxI./20):maxI)

l = legend(legStr);
set(l,'interpreter','none');