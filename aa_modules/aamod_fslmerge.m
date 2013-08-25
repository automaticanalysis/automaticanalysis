%aamod_unzipstream for taking a zipped (*.gz) file as an initial input to
%AA and turning it into a stream.
%
function [aap, resp] = aamod_fslmerge(aap, task, subjInd, sessInd)

resp='';

switch task
    case 'report'
        
    case 'doit'
    
        settings = aap.tasklist.currenttask.settings;
                
        epiFiles = aas_getfiles_bystream(aap, subjInd, sessInd, 'epi');
        
        % put files in a horizontal string to pass to fslmerge
        horizInFiles = [];
        for fileInd = 1:size(epiFiles,1)
            horizInFiles = [horizInFiles strtok(epiFiles(fileInd,:))];
        end
        
        [pth, nm, ext] = fileparts(strtok(epiFiles(1,:)));
        
        outName = fullfile(pth, ['c' nm]); % c = arbitary prefix (Combined?)
        
        cmd = sprintf('fslmerge -t %s %s', outName, horizInFiles);                

        [status, result] = aas_runfslcommand(aap, cmd);
        
        outFile = spm_select('fplist', pth, sprintf('^c%s.*', nm));
        
        aap = aas_desc_outputs(aap, subjInd, sessInd, 'epi', outFile);
        
    case 'checkrequirements'
        
        
end
