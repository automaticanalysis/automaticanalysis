%aamod_fslmaths AA module for running fslmaths commands
%
function [aap, resp] = aamod_fslmaths(aap, task, varargin)    

resp='';

switch task
    case 'report'
        
    case 'doit'
   
        settings = aap.tasklist.currenttask.settings;
        
        inStream = settings.inputstreams(1).stream;
        outStream = settings.outputstreams(1).stream;
        
        setenv('FSLOUTPUTTYPE', settings.fsloutputtype); % .nii.gz or .nii, depending on where in processing stream
        
        inFiles = cellstr(aas_getfiles_bystream(aap, varargin{:}, inStream));                
        
        outFiles = {};
         for f = 1:length(inFiles)
                thisFile = inFiles{f};
                
                [pth, nm, ext] = fileparts(thisFile);
                outFile = fullfile(pth, [settings.outputprefix nm]);
                
                cmd = sprintf(settings.cmd, thisFile, outFile);
                aas_log(aap, false, sprintf('Running: %s', cmd));
                [status, result] = aas_runfslcommand(aap, cmd);
                
                if status > 0
                    aas_log(aap, true, sprintf('Error running fslmaths: %s', result));
                end
                
                outFiles{end+1} = spm_select('fplist', fileparts(thisFile), sprintf('^%s%s.*', settings.outputprefix, nm));
         end
        
        % Describe outputs
        aap = aas_desc_outputs(aap, varargin{:}, outStream, outFiles);
        
    case 'checkrequirements'
        
        
end
