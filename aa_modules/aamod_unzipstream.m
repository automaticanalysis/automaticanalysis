%aamod_unzipstream for taking a zipped (*.gz) file as an initial input to
%AA and turning it into a stream.
%
function [aap, resp] = aamod_unzipstream(aap, task, varargin)

resp='';

switch task
    case 'report'
        
    case 'doit'
        
        settings = aap.tasklist.currenttask.settings;
        
        % Each input stream should have a corresponding output stream
        if length(settings.inputstreams) ~= length(settings.outputstreams)
            aas_log(aap, true, 'For each (zipped) input stream you need a corresponding (unzipped) output stream.');
        end
        
        if ~isfield(settings, 'split4d') || isempty(settings.split4d)
            split4d = true;        
        elseif isfield(settings, 'split4d')
            split4d = settings.split4d;        
        end
                
        % There may be more than one stream at this level, so loop through
        % each. Assume matched input and output streams (i.e. input stream
        % 1 gets unzipped to output stream 1, and so on.        
        for thisStream = 1:length(settings.inputstreams)
            
            unzippedNames = {};
            
            inStream = settings.inputstreams(thisStream).stream;
            outStream = settings.outputstreams(thisStream).stream;
            cmd = settings.unzipcmd; % typically 'gunzip -f'
            
            inFiles = cellstr(aas_getfiles_bystream(aap, varargin{:}, inStream));
            
            aas_log(aap, false, sprintf('Found %d files to unzip for stream %s.', length(inFiles), inStream));
            
            for f = 1:length(inFiles)
                thisFile = inFiles{f};
                [pth, nm, ext] = fileparts(thisFile); % ext should be .gz, so nm should be unzipped name
                unzipCmd = sprintf('%s %s', cmd, thisFile);
                system(unzipCmd);
                outFile = fullfile(pth, nm);
                
                %TODO TOFIX HACK - find any files ending in .nii-0 and
                %rename them to .nii
                if strcmp(outFile(end-5:end), '.nii-0')
                    system(sprintf('mv %s %s', outFile, outFile(1:end-2)));
                    outFile = outFile(1:end-2);
                end                                 
                
                
                % check whether we should try to split 4D to 3D                
                if split4d
                    splitFiles = aas_split4d(outFile, true); % split file and delete 4d after split
                    for thisOut = 1:size(splitFiles,1)
                        unzippedNames{end+1} = strtok(splitFiles(thisOut,:));
                    end
                else
                    unzippedNames{end+1} = outFile;
                end
                
            end
            
            % describeoutputs
            aap = aas_desc_outputs(aap, varargin{:}, outStream, unzippedNames);
            
        end % looping through streams
        
    case 'checkrequirements'
        
        
end