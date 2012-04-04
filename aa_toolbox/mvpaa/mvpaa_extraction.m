% MVPAA Extraction
% Extracts ROI data from aap.tasklist.currenttask.settingsIs

function Betas = mvpaa_extraction(aap, data, indROI, voxels)

% Check that it's worth to extract data
if sum(~isnan(data{1,1,1}(indROI))) > aap.tasklist.currenttask.settings.minVoxels
    Betas = ones(voxels, ...
        aap.tasklist.currenttask.settings.conditions, ...
        aap.tasklist.currenttask.settings.blocks, ...
        length(aap.tasklist.currenttask.settings.sessions));
    Betas = Betas*NaN;
    
    for s=1:aap.tasklist.currenttask.settings.sessions
        for b=1:aap.tasklist.currenttask.settings.blocks
            for c=1:aap.tasklist.currenttask.settings.conditions
                Betas(1:voxels,c,b,s) = data{c,b,s}(indROI);
            end
        end
    end
else
    Betas = [];
end