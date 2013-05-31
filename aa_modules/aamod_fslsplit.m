%aamod_fslmaths AA module for running fslmaths commands
%
function [aap, resp] = aamod_fslsplit(aap, task, subjInd, sessInd)    

resp='';

switch task
    case 'report'
        
    case 'doit'
           
        settings = aap.tasklist.currenttask.settings;
        
        %inStream = settings.inputstreams(1).stream;
        %outStream = settings.outputstreams(1).stream;
        
        setenv('FSLOUTPUTTYPE', settings.fsloutputtype); % .nii.gz or .nii, depending on where in processing stream
        
   
        inFiles = aas_getfiles_bystream(aap, subjInd, sessInd, 'epi');
        
        outFiles = [];
        
        for fileInd = 1:size(inFiles,1)
            thisFile = strtok(inFiles(fileInd,:));
            [pth, nm] = fileparts(thisFile);
            cmd = sprintf('fslsplit %s %s', thisFile, fullfile(pth, nm));
            [status, result] = aas_runfslcommand(aap, cmd);
            if status > 0
                aas_log(aap, true, sprintf('Error using fslsplit: %s', result));
            end
        
            outFiles = strvcat(outFiles, spm_select('fplist', pth, sprintf('^%s[0-9].*', nm)));
            
            if isempty(outFiles)
                aas_log(aap, true, sprintf('No output files found after fslsplitting %s', thisFile));
            end
        end
        

        
        % describe outputs
        aap = aas_desc_outputs(aap, subjInd, sessInd, 'epi', outFiles);
        
   case 'checkrequirements'
        
        
end