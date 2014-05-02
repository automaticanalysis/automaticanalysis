% MVPAA Extraction
% Extracts ROI data from aap.tasklist.currenttask.settingsIs

function [Betas, selcon] = mvpaa_extraction(aap, data, indROI, voxels)

% Check that it's worth to extract data
% Avoid basing decision on missing column
selcon = 1;
while sum(~isnan(data{selcon,1,1})) == 0
    selcon = selcon+1
end

if sum(~isnan(data{selcon,1,1}(indROI))) > aap.tasklist.currenttask.settings.minVoxels
    Betas = ones(voxels, ...
        aap.tasklist.currenttask.settings.conditions, ...
        aap.tasklist.currenttask.settings.blocks, ...
        aap.tasklist.currenttask.settings.sessions);
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