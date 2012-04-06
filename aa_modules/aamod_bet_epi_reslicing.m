% AA module
% Runs EPI slicing after BET
% aamod_realign should be run before running this module 

function [aap,resp]=aamod_bet(aap,task,i)

resp='';

switch task        
    case 'summary'
        subjpath=aas_getsubjpath(i);
        resp=sprintf('Align %s\n',subjpath);
        
    case 'report'
        
    case 'doit'
        
        warning off
        
        tic
        
        % RESLICE THE MASKS & MESHES
        EPIfn = aas_getfiles_bystream(aap,i,1,'meanepi');
        if size(EPIfn,1) > 1
            % Not warped!
            for a = 1:size(EPIfn,1)
                if ~strcmp(EPIfn(a,1), 'w')
                    EPIfn = EPIfn(a,:);
                    break
                end
            end
            EPIfn = EPIfn(1,:);
            fprintf('\tSeveral mean EPIs found, considering: %s\n', EPIfn)
        end
        
        fprintf('Reslicing brain masks to mean EPI\n')
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
        outMask=aas_getfiles_bystream(aap,i,'BETmask');
        spm_reslice(strvcat(EPIfn, outMask), resFlags);         
        % Get the images we resliced
        outMaskEPI = '';
        for d = 1:size(outMask,1)
            [mpth mnme mext]=fileparts(outMask(d,:));
            outMaskEPI = strvcat(outMaskEPI, fullfile(mpth,['r' mnme mext]));
        end
        
        %% DESCRIBE OUTPUTS!
        aap=aas_desc_outputs(aap,i,'epiBETmask',outMask);
        
        time_elapsed
end