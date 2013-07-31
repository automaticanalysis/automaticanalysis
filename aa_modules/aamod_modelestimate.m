function [aap, resp] = aamod_modelestimate(aap, task, subjInd)
% Estimate model that has already been specified (using, e.g.,
% aamod_firstlevel_modelspecify).

resp='';

switch task
    case 'report'
        
    case 'doit'
        
        settings = aap.tasklist.currenttask.settings;
        
        % options
        saveResids = settings.saveresids;
        
        
        % keep track of original directory so we can go back here
        startingDir = pwd;
        
        
        % Get SPM, note path, and load it
        SPMpath = aas_getfiles_bystream(aap, subjInd, 'firstlevel_spm');        
        [analysisDir, fname, ext] = fileparts(SPMpath);        
        load(SPMpath);
        SPM.swd = analysisDir;
        cd(analysisDir);
        
        
        % Estimatie the model differently depending on if we want to save
        % residuals or not        
        
        if saveResids
            
            global defaults;
            if isempty(defaults)
                spm_get_defaults;
            end
            
            % change SPM defaults to save more than 64 images
            originalMaxres = defaults.stats.maxres;
            
            defaults.stats.maxres = Inf;            
            
            spm_spm_saveresids(SPM);            
            residualImages = spm_select('fplist', analysisDir, '^ResI_.{4}\..{3}$');            
            aap = aas_desc_outputs(aap, subjInd, 'epi', residualImages);
            
            % change SPM defaults back
            defaults.stats.maxres = originalMaxres;
        else
            spm_spm(SPM);
        end
                
        % describe outputs that are common to both
        
        aap = aas_desc_outputs(aap, subjInd, 'firstlevel_spm', SPMpath);
        
        %  firstlevel_betas (includes related statistical files)
        allbetas = dir(fullfile(analysisDir,'beta_*'));
        betafns=[];
        for betaind=1:length(allbetas);
            betafns = strvcat(betafns,fullfile(analysisDir,allbetas(betaind).name));
        end
        
        otherfiles={'mask', 'ResMS', 'ResMS', 'RPV'};
        for thisFile = otherfiles
            betafns = strvcat(betafns, spm_select('fplist', analysisDir, sprintf('^%s\\..{3}$', thisFile{:})));
        end
        aap = aas_desc_outputs(aap, subjInd, 'firstlevel_betas', betafns);
        
        % add brainmask as its own
        aap = aas_desc_outputs(aap, subjInd, 'firstlevel_brainmask', spm_select('fplist', analysisDir, '^mask\..{3}$'));
        
        
        cd(startingDir);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap, 1, sprintf('Unknown task %s', task));
end