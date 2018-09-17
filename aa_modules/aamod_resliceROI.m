%aamod_resliceROI reslices ROI to match subject images
%
function [aap, resp] = aamod_resliceROI(aap, task, varargin)

resp='';

switch task
    case 'report'
        
    case 'doit'
        
        settings = aap.tasklist.currenttask.settings;
        
        % First input stream assumed to be image we are reslicing to (i.e. is
        % unchanged); following input streams are resliced and matched to the
        % first.
        %
        % There should thus be N input streams, and N-1 output streams, because all
        % of the input streams except the first will be resliced.
        
        inStreams = aap.tasklist.currenttask.inputstreams;
        outStreams = aap.tasklist.currenttask.outputstreams;
        
        
        % Get reference image from the first stream
        % try subject leve, if not, look in first session (e.g. for
        % meanepi)
        try
            refImg = aas_getfiles_bystream(aap, varargin{:}, inStreams.stream{1});
        catch
            refImg = aas_getfiles_bystream(aap, varargin{:}, 1, inStreams.stream{1});
        end
        
        
        % If more than 1, assume they are registered (e.g. EPI) and use the first
        if size(refImg,1) > 1
            refImg = strtok(refImg(1, :));
            aas_log(aap, false, sprintf('Found more than one image in first stream; using first (%s).', refImg));
        end
        
        imgs{1} = refImg;
        for streamInd = 2:length(inStreams.stream)
            imgs{streamInd} = aas_getfiles_bystream(aap, varargin{:}, inStreams.stream{streamInd});
        end
        
        
        % reslice
        flags.which = settings.which;   % generally 1, which is reslice all but first
        flags.mean = settings.mean;     % generally 0, which does not create mean image
        flags.interp = settings.interp; % generally 0 for NN, preserving binary
        flags.prefix = 'r';
        spm_reslice(strvcat(imgs), flags);
        
        
        % describe outputs (assume r prepended to input images
        for streamInd = 1:length(outStreams.stream)
            inImgs = imgs{streamInd+1}; % +1 because first img is the reference/input stream
            outImgs = {};
            for thisImg = 1:size(inImgs,1)
                inImg = strtok(inImgs(thisImg,:));
                [pth, nm, ext] = fileparts(inImg);
                outImgs{end+1} = fullfile(pth, [flags.prefix nm ext]);
            end
            
            aap = aas_desc_outputs(aap, varargin{:}, outStreams.stream{streamInd}, outImgs);
        end
        
        
        
        
    case 'checkrequirements'
        
end