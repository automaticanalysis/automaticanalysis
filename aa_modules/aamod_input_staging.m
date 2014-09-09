% aamod_input_staging
%
% A simple module that collects input streams.
%
% Useful when connecting pipelines from remote machines to miminimize data
% duplication and copying.  I have yet to think of a specific use case :) 
% Caching might do the same thing...
%
% Also can be used a choke point in the pipeline to make sure that these
% input streams are finished before spawning parallel workers that use
% them.  For example, in across-subject cross-validation we want to spawn
% an AA job for each subject, but each job requires the data from ALL 
% subjects.  Subject-domain modules don't usually wait for other subjects to
% finish, so it's possible that in a  parallel situation we are trying to
% collect data that isn't done yet.  So, having a study-domain staging
% point for that input stream ensures that the streams are up-to-date.
%
% This module was created to be a general as possible, for any combinations
% of domains and input streams.
%
% ---------------------------------------
% created by:   cwild   2014-04-20
% last updated: cwild   2014-04-20

function [aap, resp] = aamod_input_staging(aap, task, varargin)

resp='';

switch task
    case 'report'
        
    case 'doit'

        currentIndices = [varargin{:}];
        
        % Need the module index and domain
        moduleInd = aap.tasklist.currenttask.modulenumber;
        moduleDomain = aap.tasklist.currenttask.domain;
        
        % Cell arrays of stream names
        inputStreams = aap.tasklist.currenttask.inputstreams.stream;
        outputStreams = aap.tasklist.currenttask.outputstreams.stream;
              
        % Now, we pass each input stream to the corresponding output
        for iI = 1 : length(inputStreams)
            
            % First, make sure that we duplicate the input to output
           if ~ismember(inputStreams{iI}, outputStreams)
               aas_log(aap, 1, sprintf('Input stream %s not passed as output stream!', inputStreams{iI}));
           end
           
           aas_log(aap, 0, sprintf('Passing input stream %s to output.\n', inputStreams{iI}));
           
           % Domain of input stream...
           inputDomain = aap.internal.inputstreamsources{moduleInd}.stream(iI).sourcedomain;

           % Is used to get all possible locations for the stream data
           deps = aas_getdependencies_bydomain(aap, inputDomain, moduleDomain, currentIndices);
           
           % Now seach all dependency locations for the stream header 
           for dI = 1 : length(deps)
              
               % Get the path given the current dependency...
               inputPath = aas_getpath_bydomain(aap, deps{dI}{1}, deps{dI}{2});
               
               % This is where the header file should be
               inputDesc = fullfile(inputPath, sprintf('stream_%s_inputto_%s.txt', inputStreams{iI}, aap.tasklist.currenttask.name));
               
               % If exists, that means the data is there, so describe the
               % output accordingly.
               if exist(inputDesc, 'file')
                   streamFiles = aas_getfiles_bystream(aap, deps{dI}{2}, inputStreams{iI});
                   aas_desc_outputs(aap, deps{dI}{2}, inputStreams{iI}, streamFiles);
               end
               
           end
           
        end
             
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end

end

