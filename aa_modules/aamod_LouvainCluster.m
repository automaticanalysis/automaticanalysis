% aamod_LouvainClustering
%
% Details to come....

% CW - 2014-03-17

function [aap, resp] = aamod_LouvainCluster(aap, task, varargin)

resp='';

switch task
    case 'report'
        
    case 'doit'
        
        % Get the settings for this task
        settings = aap.tasklist.currenttask.settings;
        
        % By default, we're using all subjects
        subjectI = 1:length(aap.acq_details.subjects);
        
        moduleDomain = aap.tasklist.currenttask.domain;
        if ~ismember(moduleDomain, {'study' 'subject'})
            aas_log(aap, 1, sprintf('Somehow invalid module domain: %s', moduleDomain));
        end
        
        isXVal = false;
        
        indices = [varargin{:}];
        
        % Different possibilities if we're at subject level
        if strcmp(moduleDomain, 'subject')
            
            % Figure out if this a cross-validation run. We'll know because
            % the input stream will have a 'forcedomain' tag
            [modN modName modI] = aas_getmoduleindexfromtag(aap, aap.tasklist.currenttask.name);
            taskSchema = aap.schema.tasksettings.(modName)(modI);
            
            for iI = 1 : length(taskSchema.inputstreams.stream)
                if iscell(taskSchema.inputstreams.stream)  && ~isstruct(taskSchema.inputstreams.stream{1})
                    inputSchema = taskSchema.inputstreams.stream{iI};
                else
                    inputSchema = taskSchema.inputstreams.stream{1}(iI);
                end
                
                 if isstruct(inputSchema) && isfield(inputSchema.ATTRIBUTE,'forcedomain') 
                    isXVal = true;
                    break
                 end
            end

            if isXVal
                subjectI(indices) = [];   % Pop current subject if this is corss validation
            else
                subjectI = indices;       % Otherwise, we'll ONLY use the current subject
            end 
        end
        
        % The names of the output streams change depending on the .xml
        outputStreams = settings.outputstreams.stream;
        
        % Avoid issues with external libraries and matlab...
        ldpath = ['/usr/lib/x86_64-linux-gnu/libgfortran.so.3' ':/usr/lib/libhdf5.so.6'];
        setenv('LD_PRELOAD', ldpath);
        
        % Useful path variables
        pyPath = settings.toolsdir;
        binPath = fullfile(pyPath, 'bin');
        modulePath = aas_getpath_bydomain(aap, aap.tasklist.currenttask.domain, indices);
        workPath = fullfile(modulePath, 'work');
        if ~exist(workPath), mkdir(workPath); end
        
        % We have to change directories during this analysis (because of
        % how some of the bash and python scripts are coded), so let's make
        % sure we go back to where we were before.
        prevDir = pwd;
        cd(modulePath);       
        
        % Empty struct array for initialization...
        connectivityMatrices = struct('seed', {}, 'seedVoxMm', {}, 'seedVoxInd', {}, 'seedSpace', {}, 'targetNames', {}, 'targetVoxInd', {}, 'targetSpace', {}, 'correlationMatrix', {});
            
        % Collect the connectivity matrices for each subject
        aas_log(aap,false,'Loading Subjects ...');
        for subInd = subjectI
            aas_log(aap,false,sprintf('\t%s', aap.acq_details.subjects(subInd).subjname));
            matrixFile = aas_getfiles_bystream(aap, subInd, 'firstlevel_fConnMatrix_avg');
            load(matrixFile);
            
            if ~exist('avgMatrices', 'var');
                aas_log(aap, 1, sprintf('%s doesn''t contain a data structure called ''avgMatrices''', matrixFile));
            end
            
            % Could crash here if a subject has a different number of ROIs
            % than the others, but I'm too lazy to check, so instead you
            % can just read this line if it crashes here.
            connectivityMatrices(end+1, :) = avgMatrices; 
        end
        aas_log(aap,false,'');
        
        numROIs = size(connectivityMatrices, 2);
        numSubs = size(connectivityMatrices, 1);
        
        % Which ROI do we want analyze?
        if isempty(settings.whichROI)
            aas_log(aap, 1, 'settings.whichROI is empty! Which ROI are we using??');
            
        elseif isnumeric(settings.whichROI)
            
            % Could be set using an index...
            roiI = settings.whichROI;
            
            if roiI > numROIs, aas_log(aap, 1, 'There are only %d seed ROIs, but whichROI=%d', numROIs, roiI); end 
        
        else
            
            roiI = -1;
            
            % Could be specified as a string...
            for i = 1 : numROIs
               roiNames = unique({connectivityMatrices(:,i).seed});              
               m = regexp(roiNames{1}, settings.whichROI);
               if m
                   roiI = i;
                   break
               end 
            end
            
            if roiI < 0, aas_log(aap, 1, sprintf('Can''t find ROI with name %s', settings.whichROI)); end
        end
        
        % Check the ROI names, warn user if they don't correspond across subjects
        roiNames = unique({connectivityMatrices(:,roiI).seed});      
        if length(unique(roiNames)) > 1
            aas_log(aap, 0, sprintf('\nWarning: seed ROIs might be different across subjects.\nUsing these seeds: %s\n', strjoin(roiNames, ', ')), 'red');
        else
            aas_log(aap, 0, sprintf('\nWorking with seed region: %s', roiNames{1}));
        end
        
        % Trim the ROIs we aren't using
        connectivityMatrices = connectivityMatrices(:, roiI);
        
        % Check that the ROI at roiI has the same number of voxels (rows)
        % across subjects.
        numVoxPerROI = unique(arrayfun(@(x) size(x.correlationMatrix, 1), connectivityMatrices));
        if numel(numVoxPerROI) > 1
           aas_log(aap, 1, 'Seed ROIs have different number of voxels across subjects!'); 
        end
        
        % Concatenate the conectivity matrices from all subjects,
        % aligning the source voxels (rows).
        groupMatrix = [connectivityMatrices(:).correlationMatrix];
        
        % NaNs make Mark's code barf, so let's filter out NaN columns
        nanCols = any(isnan(groupMatrix));
        groupMatrix(:, nanCols) = [];
        
        aas_log(aap,false,'Performing Louvain Clustering!');
        
        louvainResults = struct('Q', [], 'Modules', []);
        
        Qs = zeros(settings.numiterations, 1); % modularity scores for each iteration
        Ms = {};                               % module memberships for each iteration
        aas_log(aap,false,sprintf('%70s', ''));
        for iter = 1 : settings.numiterations
            
            aas_log(aap,false,sprintf('%s%-70s', repmat(sprintf('\b'),1,70), sprintf('Iteration %05d / %05d', iter, settings.numiterations)));
            
            % Save the group matrix, so  subsequent scripts and exe's can read it.
            fName = sprintf('%s_%05d_groupMatrix.mat', roiNames{1}, iter);
            save(fName, 'groupMatrix');
            
            % These are the files that are written during the
            % clustering analysis.
            graphFile = fullfile(workPath, [fName '-to-eta.vtx']);
            edgeFile = fullfile(workPath, [fName '-to-eta.vtx.out.edge']);
            binFile = fullfile(workPath, [fName '-to-eta.vtx.out.bin']);
            modListFile = fullfile(workPath, [fName '-to-eta.vtx.out.coms']);
            modScoreFile = fullfile(workPath, [fName '-to-eta.vtx.out.mod']);
            
%             aas_log(aap, 0, 'Converting group connectivity matrix to VTX');
            cmd = sprintf('python %s ./%s groupMatrix CORR', fullfile(pyPath, 'mat2eta.py'), fName);
            [s w] = aas_shell(cmd);
            
%             aas_log(aap, 0, 'Computing eta^2 matrix');
            cmd = sprintf('%s %s', fullfile(binPath, 'eta'), graphFile);
            [s w] = aas_shell(cmd);
            
            % Evidently, we could run ML clustering, but Mark commented it
            % out.  Add it back in by buildign the correct system command
%             aas_log(aap, 0, 'NOT Running Machine Learning clustering')
            % #$PROD_PYTHON eta2ml.py ./work/$1-to-eta.vtx
            
%             aas_log(aap, 0, 'Converting eta^2 matrix to a graph');
            cmd = sprintf('python %s %s', fullfile(pyPath, 'eta2graph.py'), graphFile);
            [s w] = aas_shell(cmd);
            
%             aas_log(aap, 0, '\nPerforming Louvain Analysis');
            cmd = sprintf('%s -i %s -o %s', fullfile(binPath, 'convert'), edgeFile, binFile);
            [s w] = aas_shell(cmd);
            
            cmd = sprintf('%s %s -l -1 >| %s 2>| %s', fullfile(binPath, 'community'), binFile, modListFile, modScoreFile);
            [s w] = aas_shell(cmd);
            
%             aas_log(aap, 0, '\nCollecting Louvain results ...');
            cmd = sprintf('python %s %s', fullfile(pyPath, 'louvain2mat.py'), modListFile);
            [s w] = aas_shell(cmd);
            
            delete(fName);

            % Gather results
            resultsFile = [modListFile '-graph-coms.mat'];
            results = load(resultsFile);
            louvainResults(iter).Q = single(results.mod);
            louvainResults(iter).Modules = single(results.coms);
            
        end
        
        % Delete temp working files
        rmdir(workPath,'s');
        
        resultsFile = fullfile(modulePath, 'moduleResults.mat');
        save(resultsFile, 'louvainResults');
        aap = aas_desc_outputs(aap, indices, outputStreams{1}, resultsFile);
       
        aas_log(aap,false,'');
        [maxQ maxI] = max([louvainResults.Q]);
        aas_log(aap,false,sprintf('Maximum modularity obtained = %.04f\n', maxQ));
        modules = louvainResults(maxI).Modules;
        
        % Now we need info about seed voxel locations and their space, so 
        % we can create images, masks, and movies!
        voxI = connectivityMatrices(1).seedVoxInd;
        space = connectivityMatrices(1).seedSpace;

        % Store maps made from the group here
        imagePath = fullfile(modulePath, 'modulemaps');
        if ~exist(imagePath), mkdir(imagePath); end
        
        % We will be writing out individual module maps, track them here
        moduleImgs = [];
        
        % Create SPM Vol structure to use as a template for our module maps
        Vt = struct('fname', '', ...
                    'dim',   space.dim, ...
                    'mat',   space.mat, ...
                    'dt',    [4 0], ...      % Don't need much storage for a mask image
                    'pinfo', [1 0 0]', ...
                    'n',     [1 1]);

        % I think that the number of columns here is the number of
        % solutions?  Or different iterations?  Regardless, the last column
        % is the one with the fewest modules... so let's use that! :)
        numIters = size(modules, 2);

        % How many voxels in the results?
        numLabelledVox = size(modules, 1);
        
        if numLabelledVox ~= size(groupMatrix, 1)
            aas_log(aap, 1, 'Louvain results have wrong number of voxels!');
        end
        
        % Get/sort the modules by the number of voxels in each
        voxLabels = single(modules(:, numIters));
        [numVoxPerMod, modIDs] = hist(voxLabels, 0:max(voxLabels));
        [numVoxPerMod, sortI] = sort(numVoxPerMod, 'descend');
        modIDs = modIDs(sortI);
        totalNumMods = length(modIDs);
        keepModI = find(numVoxPerMod > settings.minVoxPerMod);
        numKeepMods = length(keepModI);
        
        aas_log(aap, 0, sprintf('\n%d total modules, %d of which have > %d voxels.\n', totalNumMods, numKeepMods, settings.minVoxPerMod));

        % Use this to create average "fingerprints" for each module, this
        % is different than the group matrix that we used for clustering,
        % becuase that was concatenated across subjects (here we are
        % averaging).
        avgVoxToTargetConnectivity = cat(3, connectivityMatrices.correlationMatrix);
        avgVoxToTargetConnectivity = nanmean(avgVoxToTargetConnectivity, 3);

        % Collect info about the modules
        % connectivityPattern: Module X Target matrix of average connectivity between module and
        %                      target (averaged across voxels in the module, and subjects)
        % spatialPattern: Module X voxel matrix of module membership.
        moduleInfo = struct('rank', {}, 'numVox', {}, 'connectivityPattern', {}, 'spatialPattern', {});
        
        % Go through the modules in rank order
        voxLabelsUpdated = zeros(size(voxLabels));
        for mod = 1 : totalNumMods;
            
            % Update the label to something more meaningful (i.e., rank)
            mI = find(voxLabels == modIDs(mod));
            voxLabelsUpdated(mI) = mod + 1000;
            
            % Averate across voxels in the module to get a fingerprint of connectivity
            moduleInfo(mod).connectivityPattern = nanmean(avgVoxToTargetConnectivity(mI, :), 1);
            
            % Get the spatial membership
            moduleInfo(mod).spatialPattern = zeros(1, numLabelledVox);
            moduleInfo(mod).spatialPattern(mI) = 1;
            
            moduleInfo(mod).rank = mod;
            moduleInfo(mod).numVox = length(mI);
            
            if moduleInfo(mod).numVox >= settings.minVoxPerMod
                
                % Write out an image for just this module
                Ymod = zeros(space.dim);
                Ymod(voxI(mI)) = voxLabelsUpdated(mI);
                Vmod = Vt;
                Vmod.fname = fullfile(imagePath, sprintf('M%03d_%s.nii', mod, roiNames{1}));
                Vmod.descrip = sprintf('Module ranked %03d in %s, created from group (N=%d)', mod, roiNames{1}, length(subjectI));
                spm_write_vol(Vmod, Ymod);
                moduleImgs = char(moduleImgs, Vmod.fname);
            end
        end
           
        moduleInfoFile = fullfile(modulePath, 'module_info.mat');
        save(moduleInfoFile, 'moduleInfo');
        aap = aas_desc_outputs(aap, indices, outputStreams{2}, moduleInfoFile);
        
        % Always an emptry array element first :(
        moduleImgs(1,:) = [];
        aap = aas_desc_outputs(aap, indices, outputStreams{3}, moduleImgs);
        
        % Write an image with all modules labelled
        Ymod = zeros(space.dim);
        Ymod(voxI) = voxLabelsUpdated;
        Vmod = Vt;
        Vmod.fname = fullfile(imagePath, sprintf('ALLModules_%s.nii', roiNames{1}));
        Vmod.descrip = sprintf('Labelled modules in %s, created from group (N=%d)', roiNames{1}, length(subjectI));
        spm_write_vol(Vmod, Ymod);
        aap = aas_desc_outputs(aap, indices, outputStreams{4}, Vmod.fname);               

        % Generate a movie
        if settings.diagnostic
            
            % Use the free 'colorspace' function instead of vision toolbox
            func = @(x) colorspace('RGB->Lab', x);
            colours = distinguishable_colors(numKeepMods, 'w', func);
            
            spmFig = spm_figure('FindWin', 'Graphics');
            if isempty(spmFig)
                spmFig = spm_figure('Create', 'Graphics', 'Graphics');
            end
            spm_figure('Clear','Graphics');
            
            spm_orthviews('reset');
            h = spm_orthviews('Image', aap.directory_conventions.T1template);
            for mod = 1 : numKeepMods;
                mI = find(voxLabels == modIDs(keepModI(mod)));
                [x, y, z] = ind2sub(space.dim, voxI(mI));
                spm_orthviews('AddColouredBlobs', h, [[x y z]'; ones(1, size(x,1))], ones(1, size(x,1)), space.mat, colours(mod,:), sprintf('mod%d',mod));
            end
            
            aas_checkreg_avi(aap, indices, 2, '', -80:1:50);
        end
        
        % Finally, go back to where we came from
        cd(prevDir);
        
    case 'checkrequirements'
        
        ldpath = ['/usr/lib/x86_64-linux-gnu/libgfortran.so.3' ':/usr/lib/libhdf5.so.6'];
        setenv('LD_PRELOAD', ldpath);

        [s w] = aas_shell('python -c "import numpy"', 1);
        if s ~= 0
            aas_log(aap, 1, 'Can''t load python library "numpy.py"');
        end
        
        [s w] = aas_shell('python -c "import networkx"', 1);
        if s ~= 0
            aas_log(aap, 1, 'Can''t load python library "networkx.py"');
        end
        
        [s w] = aas_shell('python -c "import h5py"', 1);
        if s ~= 0
            aas_log(aap, 1, 'Can''t load python library "h5py.py"');
        end
        
        % Check that required python files are in the tools directory
        reqdPyFiles = {'._eta2graph.py','._eta2ml.py','._louvain2mat.py','._mat2eta.py','MEG_dynamics.py','eta2graph.py','eta2ml.py','louvain2mat.py','mat2eta.py','quickset.py','thresh_methods.py','vtx.py'};
        presentFiles = dir(fullfile(aap.tasklist.currenttask.settings.toolsdir, '*.py'));
        presentFiles = {presentFiles.name};
        fileFound = ismember(reqdPyFiles, presentFiles);
        if any(~fileFound)
           aas_log(aap, 1, sprintf('Cannot find required Python scripts: %s', strjoin(reqdPyFiles(~fileFound))));
        end
        
        % Now check for other required executables (in tools/bin)
        reqdExeFiles = {'eta', 'community', 'convert', 'hierarchy'};
        presentFiles = dir(fullfile(aap.tasklist.currenttask.settings.toolsdir, 'bin', '*'));
        presentFiles = {presentFiles.name};
        fileFound = ismember(reqdExeFiles, presentFiles);
        if any(~fileFound)
           aas_log(aap, 1, sprintf('Cannot find required executables: %s', strjoin(reqdExeFiles(~fileFound))));
        end
        
    otherwise
        aas_log(aap, 1, sprintf('Unknown task %s',task));
end

end
