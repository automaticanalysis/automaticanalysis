% Create a movie

function aas_checkreg_avi(aap, p, axisDim, suffix, slicesD)

% Check if we need to make a movie (default - yes)...
if (isfield(aap.options,'diagnostic_videos') && ~aap.options.diagnostic_videos) || ... % global
        isfield(aap.tasklist.currenttask.settings,'diagnostic') && ... % module-specific
        ~isstruct(aap.tasklist.currenttask.settings.diagnostic) && ...
        ~aap.tasklist.currenttask.settings.diagnostic
    aas_log(aap,false,'INFO: Diagnostic videos disabled (check aap.options.diagnostic_videos and aap.tasksettings.<module>.diagnostic)');
else
    global st
    for v = 1:numel(st.vols)
        if isempty(st.vols{v}), break; end
        bb(:,:,v) = spm_get_bbox(st.vols{v});
    end;
    nVols = v-1;
    
    % slicesDx{1} = -85:1:85; % sagittal
    % slicesDx{2} = -120:1:90; % coronal
    % slicesDx{3} = -70:1:90; % axial
    slicesDx{1} = max(bb(1,1,:)):1:min(bb(2,1,:)); % sagittal
    slicesDx{2} = max(bb(1,2,:)):1:min(bb(2,2,:)); % coronal
    slicesDx{3} = max(bb(1,3,:)):1:min(bb(2,3,:)); % axial
    slicesInd = [];
    
    if nargin < 3
        axisDim = 2;
    end
    if nargin < 4
        suffix = '';
    end
    if nargin < 5
        if axisDim == 0 % all
            nMin = min([numel(slicesDx{1}) numel(slicesDx{2}) numel(slicesDx{3})]);
            for a = 1:3
                step = numel(slicesDx{a})/nMin;
                slicesD(a,:) = round(slicesDx{a}(1):step:slicesDx{a}(end));
            end
            slicesInd = [round(size(slicesD,2)/4) round(size(slicesD,2)/2) round(3*size(slicesD,2)/4)];
            for v = 1:nVols
                slicesImg{v} = cell(1,3);
            end;
        else slicesD = slicesDx{axisDim};
        end
    end
    
    if isempty(p)
        path = aas_getstudypath(aap);
    else
        if isnumeric(p)
            path = aas_getsubjpath(aap,p);
        else % cell
            path = aas_getpath_bydomain(aap,p{1},p{2});
        end
    end
    
    
    % Make a movie from whichever image is on SPM figure 1
    movieFilename = fullfile(path, ...
        ['diagnostic_' mfilename suffix '.avi']);
    slicesFilename0 = fullfile(path, ...
        ['diagnostic_' strrep(mfilename,'_avi','_slices') suffix]);
    
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
        figure(1);
    end
    windowSize = get(0,'ScreenSize');
    H = windowSize(4) - windowSize(2) - 100; % -100 for system menu and statusbar
    H = H - mod(H,3);
    windowSize = [1 100 H/3*2 H]; %windowSize(3) windowSize(4)];
    set(1,'Position', windowSize)
    
    for d = 1:size(slicesD,2)
        pos = [0 0 0];
        if axisDim == 0
            for a = 1:3
                pos(a) = slicesD(a,d);
            end
        else pos(axisDim) = slicesD(1,d);
        end
        spm_orthviews('reposition', pos);
        
        % Capture frame and store in aviObject
        F = getframe(1,windowSize);
        if checkmatlabreq([7;11]) % From Matlab 7.11 use VideoWriter
            writeVideo(aviObject,F);
        else
            aviObject = addframe(aviObject,F);
        end
        
        % ~FSL
        if ~isempty(find(d==slicesInd, 1))
            spm_orthviews('context_menu','Xhairs','off');
            for v = 1:nVols
                for a = 1:3
                    fr = getframe(st.vols{v}.ax{a}.ax);
                    slicesImg{v}{a} = horzcat(slicesImg{v}{a}, fr.cdata);
                end
            end
            spm_orthviews('context_menu','Xhairs','on');
        end
    end
    
    try
        % Return figure 1 to not be on top!
        fJframe.fFigureClient.getWindow.setAlwaysOnTop(false)
    catch
    end
    if checkmatlabreq([7;11]) % From Matlab 7.11 use VideoWriter
        close(aviObject);
    else
        junk = close(aviObject);
    end
    
    % ~FSL
    if exist('slicesImg','var')
        for v = 1:numel(slicesImg)
            slicesFilename = sprintf('%s_%d.jpg',slicesFilename0,v);
            img = slicesImg{v}{3};
            img(1:size(slicesImg{v}{2},1),end+1:end+size(slicesImg{v}{2},2),:) = slicesImg{v}{2};
            img(1:size(slicesImg{v}{1},1),end+1:end+size(slicesImg{v}{1},2),:) = slicesImg{v}{1};
            f = figure;
            set(f,'Position',[1 1 size(img,2) size(img,1)],'PaperPositionMode','auto','InvertHardCopy','off');
            try
                imshow(img,'Border','tight');
                set(f,'Renderer','zbuffer');
                print(f,'-djpeg','-r150',slicesFilename);
            catch
            end;
            close(f);
        end
    end
end
end
