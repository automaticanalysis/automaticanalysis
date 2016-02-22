% AA module - perform MarsBaR ROI analysis on estimated SPM design

function [aap, resp] = aamod_marsbar(aap, task, subjInd)

resp = '';

switch task
    case 'report'
        
    case 'doit'

        if ~marsbar('is_started'), marsbar; end
        
        settings = aap.tasklist.currenttask.settings;
        
        switch settings.roidef
            case 'label'
                settings.roidef = 'i';
            case 'cluster'
                settings.roidef = 'c';
            otherwise
                aas_log(aap, 1, sprintf('Invalid ROI type setting! ''%s''', settings.roitype));
        end
        

        % Initialize the output data structures
        marsbarResults = {}; % will hold contrast values, stats, etc.
        pctSignal = struct('roi', {}, 'rows', {}, 'columns',{}, 'eventSpecs', {}, 'PSC', []); % will hold percent signal change for all events in the model  

        % Clear the marsbar directory, because we generate ROIs then scan 
        % to pick them up. We don't want to include out-of-date ROIs.
        subjDir = aas_getsubjpath(aap, subjInd);
        mbDir = fullfile(subjDir, 'MarsBaR');   
        if exist(mbDir), rmdir(mbDir,'s'); end
        mkdir(mbDir);
        
        % load the SPM file as a MarsBaR design object
        spmFile = aas_getfiles_bystream(aap, subjInd, 'firstlevel_spm');
        mbSPM = mardo(spmFile);
        
        % also load the estimated SPM file, need for contrasts
        SPM = des_struct(mbSPM);
        
        if ~isfield(SPM, 'xCon') || isempty(SPM.xCon)
           aas_log(aap, 1, 'SPM model contains no contrasts! Try adding aamod_firstlevel_contrasts to your tasklist before this module.'); 
        end
        
        % Get roi imgs
        roiImgs = cellstr(aas_getfiles_bystream(aap, subjInd, 'rois'));
        allROIs = {};
        
        for r = 1 : size(roiImgs, 1)
            
            [roiPath, roiImg, roiExt] = fileparts(roiImgs{r});
            
            if settings.binarizeROIs
                fprintf('Binarizing ROI %s\n', roiImg);
                v = spm_vol(roiImgs{r});
                y = spm_read_vols(v);
                y(y~=0) = 1;
                v.fname = fullfile(roiPath, ['b' roiImg roiExt]);
                spm_write_vol(v, y);
                roiImgs{r} = v.fname;
            end
            
            % Convert image to MarsBar type ROI
            mars_img2rois(roiImgs{r}, mbDir, roiImg, settings.roidef);
            roiFiles = dir(fullfile(mbDir, sprintf('%s*.mat', roiImg)));
            fprintf('\nMarsBaR created %d ROIs from %s.\n\n', numel(roiFiles), roiImg);   
            roiFiles = fullfile(mbDir, {roiFiles.name});
            allROIs = [allROIs roiFiles];
            rois = maroi(roiFiles);
            
            % Extract the data from the ROI
            Y = aa_get_marsy(rois, mbSPM, settings.summaryfn);
            
            % MarsBaR estimation
            mbSPM = estimate(mbSPM, Y);
            mbSPM = add_contrasts(mbSPM, SPM.xCon);
            marsbarResults{r} = compute_contrasts(mbSPM);
            
            % Event definitions from the model
            [eSpecs, eNames] = event_specs(mbSPM);
            numEvents = size(eSpecs, 2);

            roiLabels = {};
            for i = 1 : numel(rois)
                roiLabels{i} = label(rois{i});
            end
            
            % init struct
            pctSignal(r) = struct('roi', roiImg, 'rows', {roiLabels}, 'columns', {eNames}, 'eventSpecs', eSpecs, 'PSC', []);
            
            % Get the percent signal esimate for all events
            for e = 1 : numEvents
                
                % Get the average duration of this event (convert to secs)
                [junk, durations] = event_onsets(mbSPM, eSpecs(:, e));
                avgDur = mean(durations .* SPM.xY.RT);
                
                % Calculate % Signal change for this event
                pctSignal(r).PSC(:,e) = event_signal(mbSPM, eSpecs(:, e), avgDur);
            end
            
        end
        
%         roiFiles(1,:) = []; % pop the first empty row of the char array 
        aap = aas_desc_outputs(aap, subjInd, 'marsbar_rois', char(allROIs));
        
        dataFile = fullfile(mbDir, 'marsbar_results.mat');
        save(dataFile, 'marsbarResults', 'pctSignal');
        aap = aas_desc_outputs(aap, subjInd, 'marsbar_results', dataFile);
        
    case 'checkrequirements'
        
        r = which('marsbar');
        if isempty(r)
            aas_log(aap, 1, 'MarsBaR is not present in the MATLAB path, please add it in your user script.');
        end
end

end

% MarsBaR uses a horrible method of handling multiple input ROIs! Rather
% than having one input that contains multiple ROIs, MarBaR has a variable
% number of input arguments, so you are supposed to specify ROIs like
% this: y = get_marsy(roi1, roi2, roi3, ..., SPM, summary);
%
% Normally, you would build a string for that command execute it with
% eval(), but we are avoiding eval() in aa because you can't compile matlab
% code that contains that function. My solution for now is to create
% overload aa_get_marsy() which has a switch statemebt for as many ROIs 
% as we think we need. OBviously this is an impractical solution if you 
% have many ROIs in the same image file.

function Y = aa_get_marsy(rois, SPM, summaryfn)
switch numel(rois)
    case 1
        Y = get_marsy(rois{1}, SPM, summaryfn);
    case 2
        Y = get_marsy(rois{1}, rois{2}, SPM, summaryfn);
    case 3
        Y = get_marsy(rois{1}, rois{2}, rois{3}, SPM, summaryfn);
    case 4
        Y = get_marsy(rois{1}, rois{2}, rois{3}, rois{4}, SPM, summaryfn);
    case 5
        Y = get_marsy(rois{1}, rois{2}, rois{3}, rois{4}, rois{5}, SPM, summaryfn);
    case 6
        Y = get_marsy(rois{1}, rois{2}, rois{3}, rois{4}, rois{5}, rois{6}, SPM, summaryfn);
    case 7
        Y = get_marsy(rois{1}, rois{2}, rois{3}, rois{4}, rois{5}, rois{6}, rois{7}, SPM, summaryfn);
    case 8
        Y = get_marsy(rois{1}, rois{2}, rois{3}, rois{4}, rois{5}, rois{6}, rois{7}, rois{8}, SPM, summaryfn);
    case 9
        Y = get_marsy(rois{1}, rois{2}, rois{3}, rois{4}, rois{5}, rois{6}, rois{7}, rois{8}, rois{9}, SPM, summaryfn);
    case 10
        Y = get_marsy(rois{1}, rois{2}, rois{3}, rois{4}, rois{5}, rois{6}, rois{7}, rois{8}, rois{9}, rois{10}, SPM, summaryfn);
    case 11
        Y = get_marsy(rois{1}, rois{2}, rois{3}, rois{4}, rois{5}, rois{6}, rois{7}, rois{8}, rois{9}, rois{10}, rois{11}, SPM, summaryfn);
    case 12
        Y = get_marsy(rois{1}, rois{2}, rois{3}, rois{4}, rois{5}, rois{6}, rois{7}, rois{8}, rois{9}, rois{10}, rois{11}, rois{12}, SPM, summaryfn);
    case 13
        Y = get_marsy(rois{1}, rois{2}, rois{3}, rois{4}, rois{5}, rois{6}, rois{7}, rois{8}, rois{9}, rois{10}, rois{11}, rois{12}, rois{13}, SPM, summaryfn);
    case 14
        Y = get_marsy(rois{1}, rois{2}, rois{3}, rois{4}, rois{5}, rois{6}, rois{7}, rois{8}, rois{9}, rois{10}, rois{11}, rois{12}, rois{13}, rois{14}, SPM, summaryfn);
    case 15
        Y = get_marsy(rois{1}, rois{2}, rois{3}, rois{4}, rois{5}, rois{6}, rois{7}, rois{8}, rois{9}, rois{10}, rois{11}, rois{12}, rois{13}, rois{14}, rois{15}, SPM, summaryfn);
    case 16
        Y = get_marsy(rois{1}, rois{2}, rois{3}, rois{4}, rois{5}, rois{6}, rois{7}, rois{8}, rois{9}, rois{10}, rois{11}, rois{12}, rois{13}, rois{14}, rois{15}, rois{16}, SPM, summaryfn);
    case 17
        Y = get_marsy(rois{1}, rois{2}, rois{3}, rois{4}, rois{5}, rois{6}, rois{7}, rois{8}, rois{9}, rois{10}, rois{11}, rois{12}, rois{13}, rois{14}, rois{15}, rois{16}, rois{17}, SPM, summaryfn);
    case 18
        Y = get_marsy(rois{1}, rois{2}, rois{3}, rois{4}, rois{5}, rois{6}, rois{7}, rois{8}, rois{9}, rois{10}, rois{11}, rois{12}, rois{13}, rois{14}, rois{15}, rois{16}, rois{17}, rois{18}, SPM, summaryfn);
    case 19
        Y = get_marsy(rois{1}, rois{2}, rois{3}, rois{4}, rois{5}, rois{6}, rois{7}, rois{8}, rois{9}, rois{10}, rois{11}, rois{12}, rois{13}, rois{14}, rois{15}, rois{16}, rois{17}, rois{18}, rois{19}, SPM, summaryfn);
    case 20
        Y = get_marsy(rois{1}, rois{2}, rois{3}, rois{4}, rois{5}, rois{6}, rois{7}, rois{8}, rois{9}, rois{10}, rois{11}, rois{12}, rois{13}, rois{14}, rois{15}, rois{16}, rois{17}, rois{18}, rois{19}, rois{20}, SPM, summaryfn);
    otherwise
        aas_log(aap, 1, 'Unsupported number of ROIs, try adding a new case to this function');
end
end

% % Use this to generate all the switch cases:
% s1 = 'rois{%d}, ';
% for i = 1 : 20
%     fprintf('case %d\n', i);
%     fprintf('    Y = get_marsy(');
%     fprintf(s1, [1:i]);
%     fprintf('SPM, summaryfn);\n');
% end


