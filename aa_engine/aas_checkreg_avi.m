% Create a movie

function aas_checkreg_avi(aap, p, axisDim, suffix, slicesD)


% Check if we need to make a movie...
if aap.tasklist.currenttask.settings.diagnostic
     warning('OFF', 'MATLAB:getframe:RequestedRectangleExceedsFigureBounds')
    
    if nargin < 3
        axisDim = 2;
    end
    if nargin < 4
        suffix = '';
    end
    if nargin < 5
        if axisDim == 1
            slicesD = -85:1:85;
        elseif axisDim == 2
            slicesD = -120:1:90;
        elseif axisDim == 3
            slicesD = -70:1:90;
        end
        
    end
    
    if isempty(p)
        path = aas_getstudypath(aap);
    else
        path = aas_getsubjpath(aap,p);
    end
    
    
    % Make a movie from whichever image is on SPM figure 1
    movieFilename = fullfile(path, ...
        ['diagnostic_' strrep(mfilename,'_avi','') suffix '.avi']);
    
    % Create movie file by defining aviObject
    if exist(movieFilename,'file')
        delete(movieFilename);
    end
    if checkmatlabreq([7;11]) % From Matlab 7.11 use VideoWriter
        aviObject = VideoWriter(movieFilename);
        open(aviObject);
    else
        aviObject = avifile(movieFilename,'compression','none');
    end
    try
        % Try to set figure 1 to be on top!
        fJframe = get(1, 'JavaFrame');
        fJframe.fFigureClient.getWindow.setAlwaysOnTop(true)
    catch
    end
    %windowSize = get(1,'Position');
    % It does not work if it's larger than the window, conservative...
    windowSize = [1 1 400 600]; %windowSize(3) windowSize(4)];
    set(1,'Position', windowSize)
    
    for d = slicesD
        if axisDim == 1
            spm_orthviews('reposition', [d 0 0])
        elseif axisDim == 2
            spm_orthviews('reposition', [0 d 0])
        elseif axisDim == 3
            spm_orthviews('reposition', [0 0 d])
        end
        
        % Capture frame and store in aviObject
        F = getframe(1,windowSize);
        %         figure(5); imagesc(F.cdata) % DEBUG CODE TO SEE IMAGE...
        if checkmatlabreq([7;11]) % From Matlab 7.11 use VideoWriter
            writeVideo(aviObject,F);
        else
            aviObject = addframe(aviObject,F);
        end
    end
    
    try
        % Return figure 1 to not be on top!
        fJframe.fFigureClient.getWindow.setAlwaysOnTop(false)
    catch
    end
    close(aviObject);
     warning('ON', 'MATLAB:getframe:RequestedRectangleExceedsFigureBounds')
end
end