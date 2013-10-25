% Create a movie
% aas_image_avi(imageFN, outlineFN, movieFN, axisDim, frameSize, rotations, outlineType)
% If movieFN is empty, then the video will be shown, not saved...
% Options for outline type:
% 'canny', 'sobel', 'prewitt', 'roberts', 'log', 'zerocross' = Outlines the volumes gently
% 'none' = No outline, but fills with small points, so you can see the
% background (good for checking segmentations...)
% 'fill' = No outline, but fills with large squares, so you cannot easily
% see background colour...
function aas_image_avi(imageFN, outlineFN, movieFN, axisDim, frameSize, rotations, outlineType)

if ischar(imageFN)
    imageFN = strvcat2cell(imageFN);
end
if nargin < 2
    outlineFN = [];
elseif ischar(outlineFN)
    outlineFN = strvcat2cell(outlineFN);
end
if nargin < 3
    [junk, movieFN] = fileparts(imageFN{1});
    
    % Make a movie file
    movieFN = fullfile(getHome, ...
        [movieFN '.avi']);
end
if nargin < 4 || isempty(axisDim)
    axisDim = 1;
end
if nargin < 5 || isempty(frameSize)
    frameSize = [400 500];
end
if nargin < 6 || isempty(rotations)
    rotations = 0;
end
if nargin < 7 || isempty(outlineType)
    outlineType = 'canny';
    % Can be 'none' (no outline...)
    % or 'sobel', 'prewitt', 'roberts', 'log', 'zerocross', 'canny'
end

% Create movie file by defining aviObject
if ~isempty(movieFN)
    if exist(movieFN,'file')
        delete(movieFN);
    end
    
    aviObject = avifile(movieFN,'compression','none');
end

% Get the figure!
figure(2)
try
    % Try to set figure 1 to be on top!
    fJframe = get(2, 'JavaFrame');
    fJframe.fFigureClient.getWindow.setAlwaysOnTop(true)
catch
end

% This does not work if it's larger than the window, conservative...
windowSize = [1 1 frameSize(1) frameSize(2)];
set(2,'Position', windowSize)

Y = cell(size(imageFN));
% Load the image
for f = 1:length(imageFN)
    Y{f} = spm_read_vols(spm_vol(imageFN{f}));
    limsY{f} = [min(Y{f}(:)) max(Y{f}(:))];
end

% Load the outline
if ~isempty(outlineFN)
    colorsB = {'r' 'g' 'b' 'c' 'm' 'y' 'w'};
    % Variables we need for outlining...
    thresh = cell(size(outlineFN));
    outlineSlice = cell(size(outlineFN));
    max_oY = -Inf;
    oY = cell(size(outlineFN));
    
    for o = 1:length(outlineFN)
        oY{o} = spm_read_vols(spm_vol(outlineFN{o}));
        if ~strcmp(outlineType, 'none') && ~strcmp(outlineType, 'fill')
            % Get a good threshold
            thresh{o} = zeros(size(Y{1},axisDim),2);
            for d = 1:size(Y{1},axisDim)
                if axisDim == 1
                    outlineSlice{o} = squeeze(oY{o}(d,:,:));
                elseif axisDim == 2
                    outlineSlice{o} = squeeze(oY{o}(:,d,:));
                elseif axisDim == 3
                    outlineSlice{o} = squeeze(oY{o}(:,:,d));
                end
                
                [outlineSlice{o} thresh{o}(d,:)] = edge(outlineSlice{o}, outlineType);
            end
            thresh{o} = mean(thresh{o});
        else
            max_oY = max(max_oY, oY{o});
        end
    end
end
% Ensure that we don't get overlap between our masks...
if strcmp(outlineType, 'none') || strcmp(outlineType, 'fill')
    for o = 1:length(outlineFN)
        oY{o}(oY{o}<max_oY) = 0;
    end
end

colormap gray

for d = 1:size(Y{1},axisDim)
    for f = 1:length(imageFN)
        h = subplot(1, length(imageFN), f);
        
        % Get image slice to draw
        if axisDim == 1
            imageSlice = squeeze(Y{f}(d,:,:));
        elseif axisDim == 2
            imageSlice = squeeze(Y{f}(:,d,:));
        elseif axisDim == 3
            imageSlice = squeeze(Y{f}(:,:,d));
        end
        
        % If present, get outline slice to draw
        if ~isempty(outlineFN)
            for o = 1:length(outlineFN)
                if axisDim == 1
                    outlineSlice{o} = squeeze(oY{o}(d,:,:));
                elseif axisDim == 2
                    outlineSlice{o} = squeeze(oY{o}(:,d,:));
                elseif axisDim == 3
                    outlineSlice{o} = squeeze(oY{o}(:,:,d));
                end
                if ~strcmp(outlineType, 'none') && ~strcmp(outlineType, 'fill')
                    outlineSlice{o} = edge(outlineSlice{o}, outlineType, thresh{o});
                end
            end
        end
        
        % Rotate slices
        imageSlice = rot90(imageSlice, rotations);
        if ~isempty(outlineFN)
            for o = 1:length(outlineFN)
                outlineSlice{o} = rot90(outlineSlice{o},rotations - 1);
            end
        end
        
        % Draw slices
        imagescnan(imageSlice, 'NanColor', [1 1 1])
        if ~isempty(outlineFN)
            hold on
            for o = 1:length(outlineFN)
                [x y] = find(flipdim(outlineSlice{o},2));
                if strcmp(outlineType, 'none')
                    scatter(x,y,round(min(frameSize./size(outlineSlice{o}))),colorsB{o}, '.')
                elseif strcmp(outlineType, 'fill')
                    scatter(x,y,round(min(frameSize./size(outlineSlice{o}))),colorsB{o}, 's')
                else
                    scatter(x,y,round(min(frameSize./size(outlineSlice{o}))),colorsB{o}, 'd')
                end
            end
            hold off
        end
        
        caxis([limsY{f}(1), limsY{f}(2)])
        axis equal off
        zoomSubplot(h, 1.2)
    end
    pause(0.01)
    drawnow
    
    if ~isempty(movieFN)
        % Capture frame and store in aviObject
        F = getframe(2,windowSize);
        aviObject = addframe(aviObject,F);
    end
end
% Save video
if ~isempty(movieFN)
    aviObject = close(aviObject);
end
try
    % Return figure 1 to not be on top!
    fJframe.fFigureClient.getWindow.setAlwaysOnTop(false)
catch
end