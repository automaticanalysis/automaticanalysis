% AA module
% Runs BET (FSL Brain Extration Toolbox) on structural
% [For correct functionality, it is recommended you run this after
% realignment and before writing the normalised EPI image
% If you do it before estimating the normalisation, make sure you normalise
% to a scull-stripped template, if at all possible!]

function [aap,resp]=aamod_bet(aap,task,p)

resp='';

switch task
    case 'domain'
        resp='subject';  % this module needs to be run once per subject
        
    case 'description'
        resp='SPM5 align';
        
    case 'summary'
        subjpath=aas_getsubjpath(p);
        resp=sprintf('Align %s\n',subjpath);
        
    case 'report'
        
    case 'doit'
        
        tic
        
        % Let us use the native space...
        Sfn = aas_getfiles_bystream(aap,p,'structural');
        EPIfn = aas_getfiles_bystream(aap,p,1,'meanepi');
        
        % Cheap and cheerful way of ensuring only one file is considered!
        if size(Sfn,1) > 1
            for a = 1:size(Sfn,1)
                % Not warped or betted!
                if ~strcmp(Sfn(a,1), 'w') && ~strcmp(Sfn(a,1), 'b')
                    Sfn = Sfn(a,:);
                    break
                end
            end
            fprintf('\tSeveral structurals found, considering: %s\n', Sfn)
        end
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
        
        [pth nme ext]=fileparts(Sfn);
        
        outStruct=fullfile(pth,['bet_' nme ext]);
        % Run BET [-R Using robust setting to avoid neck!]   
        fprintf('Run BET\n')
        [~, w]=aas_runfslcommand(aap, ...
            sprintf('bet %s %s -f %f -v',Sfn,outStruct, ...
            aap.tasklist.currenttask.settings.bet_f_parameter));
        
        % This outputs last radius from recursive command...
        indxS = strfind(w, 'radius');
        indxS = indxS(end) + 7;
        indxE = strfind(w(indxS:end), ' mm');
        indxE = indxE(1) - 2;
        SRad = w(indxS:indxS+indxE);
        
        % We don't extract the centre of gravity from here, since it needs
        % to be input in voxels... Instead get it from betted image
        Y = spm_read_vols(spm_vol(outStruct));
        Y = Y > 0;
        indY = find(Y);
        [subY_x subY_y subY_z] = ind2sub(size(Y), indY);
        COG = [mean(subY_x), mean(subY_y), mean(subY_z)];
        
        fprintf('\t...calculated c-o-g (vox): %0.4f %0.4f %0.4f  and radius (mm): %s\n', ...
            COG(1), COG(2), COG(3), SRad)
        
        %% BRAIN MASK (slightly different from inskull)
        
        V = spm_vol(fullfile(pth,['bet_' nme ext]));
        mY = spm_read_vols(V);
        
        % Mask out non-brain
        mY = mY > 0;
        
        % Then write out actual BET mask
        V.fname = fullfile(pth, ['bet_' nme '_brain_mask' ext]);
        spm_write_vol(V,mY);
        
        %% RESLICE THE MASKS & MESHES
        
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
        
        % Get the images to reslice
        D = dir(fullfile(pth, 'bet*mask*'));
        outMask = '';
        for d = 1:length(D)
            outMask = strvcat(outMask, fullfile(pth, D(d).name));
        end
        
        % Reslice
        spm_reslice(strvcat(EPIfn, outMask), resFlags); 
        
        % Get the images we resliced
        D = dir(fullfile(pth, 'rbet*mask*'));
        outMaskEPI = '';
        for d = 1:length(D)
            outMaskEPI = strvcat(outMaskEPI, fullfile(pth, D(d).name));
        end
        
        %% DIAGNOSTIC IMAGE
        % Save graphical output to common diagnostics directory
        if ~exist(fullfile(aap.acq_details.root, 'diagnostics'), 'dir')
            mkdir(fullfile(aap.acq_details.root, 'diagnostics'))
        end
        [~, mriname] = fileparts(aas_getsubjpath(aap,p));
        try
            %% Draw structural image...
            spm_check_registration(Sfn)
            
            indx = 0;
            
            % Colour the brain extracted bit pink
            spm_orthviews('addcolouredimage',1,outStruct, [0.9 0.4 0.4])
            spm_orthviews('reposition', [0 0 0])
            
            figure(spm_figure('FindWin'));
            print('-djpeg','-r75',fullfile(aap.acq_details.root, 'diagnostics', ...
                [mfilename '__' mriname '.jpeg']));
        catch
        end
        
        %% DESCRIBE OUTPUTS!
        
        % Structural image after BETting
        aap=aas_desc_outputs(aap,p,'structural', strvcat(... 
            aas_getfiles_bystream(aap,p,'structural'), ...
            outStruct));
        aap=aas_desc_outputs(aap,p,'BETmask',outMask);
        aap=aas_desc_outputs(aap,p,'epiBETmask',outMaskEPI);
        
        time_elapsed
end