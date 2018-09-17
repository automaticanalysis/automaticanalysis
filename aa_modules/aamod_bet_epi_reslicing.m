% AA module
% Runs EPI slicing after BET
% aamod_realign should be run before running this module 

function [aap,resp]=aamod_bet_epi_reslicing(aap,task,subj)

resp='';

switch task        
    case 'summary'
        subjpath=aas_getsubjpath(subj);
        resp=sprintf('Align %s\n',subjpath);
        
    case 'report'
        
    case 'doit'
        
        warning off
        
        % RESLICE THE MASKS & MESHES
        EPIimg = aas_getfiles_bystream(aap,subj,'meanepi');
        
        % With the mean EPI, we just use the first one (there really should be only one)
        if size(EPIimg,1) > 1
            EPIimg = deblank(EPIimg(1,:));
            aas_log(aap,false,sprintf('\tWARNING: Several mean EPIs found, considering: %s', EPIimg))
        end
        
        aas_log(aap,false,'Reslicing brain masks to mean EPI')
        % Get realignment defaults
        defs = aap.spm.defaults.realign;

        % Flags to pass to routine to create resliced images
        % (spm_reslice)
        resFlags = struct(...
            'interp', defs.write.interp,...       % interpolation type
            'wrap', defs.write.wrap,...           % wrapping info (ignore...)
            'mask', defs.write.mask,...           % masking (see spm_reslice)
            'which', 1,...     % what images to reslice
            'mean', 0);           % write mean image
   
        % Get files to reslice
        outMask=aas_getfiles_bystream(aap,subj,'BETmask');
        
        spm_reslice(strvcat(EPIimg, outMask), resFlags);
        
        % Get the images we resliced
        outMaskEPI = '';
        for d = 1:size(outMask,1)
            [mpth mnme mext]=fileparts(outMask(d,:));
            outMaskEPI = strvcat(outMaskEPI, fullfile(mpth,['r' mnme mext]));
            % Additionally, convert into a binary image...
            img2mask(fullfile(mpth,['r' mnme mext]))
        end
        
        %% DESCRIBE OUTPUTS!
        aap=aas_desc_outputs(aap,subj,'epiBETmask',outMaskEPI);
end
