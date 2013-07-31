%aamod_coreg_structural2template coregister structural to template
%
% This may help segmentation by improving starting estimates. Optionally,
% bring along other images (e.g. coregistered EPIs) for the ride.
%
% As with all coregistration and normalization, it's a good idea to
% manually check the results to make sure it worked!
%
function [aap, resp] = aamod_coreg_structural2template(aap, task, subjInd)

resp='';

switch task
    case 'report'
        
    case 'doit'
        
        settings = aap.tasklist.currenttask.settings;        
        inStreams = aap.tasklist.currenttask.inputstreams; % assume 1st is T1!
        outStreams = aap.tasklist.currenttask.outputstreams;
        
        
        if isempty(settings.template)
            template = fullfile(spm('Dir'), 'canonical', 'avg152T1.nii');
        end
        
        if ~exist(template, 'file')
            aas_log(aap, true, sprintf('Specified template %s not found.', template));
        end
        
        Vref = spm_vol(template);
        
        imgs = {};
            
        % Get files to update - little bit of a hack to get session-wise
        % for EPI, even though this is run at the subject level.
        for streamInd = 1:length(inStreams.stream)
            thisStream = inStreams.stream{streamInd};
            
            if strcmp(thisStream, 'epi')                
                for sess=aap.acq_details.selected_sessions
                    imgs{streamInd} = [];
                    imgs{streamInd} = strvcat(imgs{streamInd}, aas_getfiles_bystream(aap, subjInd, sess, thisStream));
                end
            else                
                imgs{streamInd} = aas_getfiles_bystream(aap, subjInd, thisStream);
            end
        end
        
        Vest = spm_vol(strtok(imgs{1})); % estimate based on first image
        
        x = spm_coreg(Vref, Vest); % TODO: add flags option? 
        M  = inv(spm_matrix(x));
        
        
        % now apply transformations, and describe outputs (by session for
        % EPI)
        for streamInd = 1:length(outStreams.stream)
            thisStream = outStreams.stream{streamInd};
            
            if ismember(thisStream, {'epi'})                
                for sess=aap.acq_details.selected_sessions
                    imgs = aas_getfiles_bystream(aap, subjInd, sess, thisStream);
                    
                    for imgInd = 1:size(imgs,1)
                        thisImg = strtok(imgs(imgInd,:));
                        spm_get_space(thisImg, M*spm_get_space(thisImg));
                    end
                    
                    aap = aas_desc_outputs(aap, subjInd, sess, thisStream, imgs);
                end % going through sessions
                
            else
                imgs = aas_getfiles_bystream(aap, subjInd, thisStream);
                
                for imgInd = 1:size(imgs,1)
                    thisImg = strtok(imgs(imgInd,:));
                    spm_get_space(thisImg, M*spm_get_space(thisImg));
                end
                
                aap = aas_desc_outputs(aap, subjInd, thisStream, imgs);
            end
            
            
        end
    case 'checkrequirements'
        
end