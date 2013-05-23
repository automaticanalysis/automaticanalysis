%aamod_unzipstream for taking a zipped (*.gz) file as an initial input to
%AA and turning it into a stream.
%
function [aap, resp] = aamod_unzipstream(aap, task, varargin)


resp='';

switch task
    case 'report'
    case 'doit'
        
        % There may be more than one stream at this level, so loop through
        % each. Assume matched input and output streams (i.e. input stream
        % 1 gets unzipped to output stream 1, and so on.
        
        for thisStream = 1:length(aap.tasklist.currenttask.settings.inputstreams)
            
            unzippedNames = {};
            
            inStream = aap.tasklist.currenttask.settings.inputstreams(thisStream).stream;
            outStream = aap.tasklist.currenttask.settings.outputstreams(thisStream).stream;
            cmd = aap.tasklist.currenttask.settings.unzipcmd; % typically 'gunzip'
            
            inFiles = cellstr(aas_getfiles_bystream(aap, varargin{:}, inStream));
            
            aas_log(aap, false, sprintf('Found %d files to unzip for stream %s.', length(inFiles), inStream));
            
            for f = 1:length(inFiles)
                thisFile = inFiles{f};
                [pth, nm, ext] = fileparts(thisFile); % ext should be .gz, so nm should be unzipped name
                unzipCmd = sprintf('%s %s', cmd, thisFile);
                system(unzipCmd);
                unzippedNames{end+1} = fullfile(pth, nm);
            end
            
            % describeoutputs
            aap = aas_desc_outputs(aap, varargin{:}, outStream, unzippedNames);            
        end % looping through streams
        
    case 'checkrequirements'
        
        
end