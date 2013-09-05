% AA module
% Runs BET (FSL Brain Extration Toolbox) on structural (usually)
% [For best functionality, it is recommended you run this after
% realignment and before writing the normalised EPI image
% If you do it before estimating the normalisation, make sure you normalise
% to a scull-stripped template, if at all possible!]

function [aap,resp]=aamod_bet(aap,task,subj)

resp='';

switch task
    case 'report'
        
    case 'doit'
        
        % Find out what stream we should BET
        inputstream = aap.tasklist.currenttask.inputstreams.stream;
        % And the names of the output streams
        outputstream = aap.tasklist.currenttask.outputstreams.stream;
        % And which are the streams which we output...
        outputstream = outputstream(~[strcmp(inputstream,outputstream)]);
        
        % Get the image...
        Simg = aas_getfiles_bystream(aap,subj,inputstream{:});
        
        % Which file is considered, as determined by the structural parameter!
        if size(Simg,1) > 1
            Simg = deblank(Simg(aap.tasklist.currenttask.settings.structural, :));
            fprintf('WARNING: Several %s found, considering: \n', inputstream{:})
            for t = 1:length(aap.tasklist.currenttask.settings.structural)
                fprintf('\t%s\n', Simg(t,:))
            end
        end
        
        % Image that we will be using for BET...
        cSimg = deblank(Simg(1,:));
        [Spth Sfn Sext]=fileparts(cSimg);
        
        % Structural image (after BETting, if we mask...)
        bSimg=fullfile(Spth,['bet_' Sfn Sext]);
        
        if aap.tasklist.currenttask.settings.robust
            % Run BET [-R Using robust setting to improve performance!]
            fprintf('1st BET pass (recursive) to find optimal centre of gravity and radius\n')
            [junk, w]=aas_runfslcommand(aap, ...
                sprintf('bet %s %s -f %f -v -R',cSimg,bSimg, ...
                aap.tasklist.currenttask.settings.bet_f_parameter));
        else
            fprintf('1st BET pass\n')
            [junk, w]=aas_runfslcommand(aap, ...
                sprintf('bet %s %s -f %f -v ',cSimg,bSimg, ...
                aap.tasklist.currenttask.settings.bet_f_parameter));
        end
        
        % This outputs last radius from recursive command...
        indxS = strfind(w, 'radius');
        indxS = indxS(end) + 7;
        indxE = strfind(w(indxS:end), ' mm');
        indxE = indxE(1) - 2;
        SRad = w(indxS:indxS+indxE);
        
        % We don't extract the centre of gravity from here, since it needs
        % to be input in voxels... Instead get it from betted image
        Y = spm_read_vols(spm_vol(bSimg));
        Y = Y > 0;
        indY = find(Y);
        [subY_x subY_y subY_z] = ind2sub(size(Y), indY);
        COG = [mean(subY_x), mean(subY_y), mean(subY_z)];
        
        fprintf('\t...calculated c-o-g (vox): %0.4f %0.4f %0.4f  and radius (mm): %s\n', ...
            COG(1), COG(2), COG(3), SRad)
        
        if aap.tasklist.currenttask.settings.masks
            fprintf('2nd BET pass extracting brain masks \n')
            % Run BET [-A Now let's get the brain masks and meshes!!]
            [junk, w]=aas_runfslcommand(aap, ...
                sprintf('bet %s %s -f %f -c %0.4f %0.4f %0.4f -r %s -v -A',cSimg,bSimg, ...
                aap.tasklist.currenttask.settings.bet_f_parameter, COG(1), COG(2), COG(3), SRad)...
                );
            
            % This outputs last radius from recursive command...
            indxS = strfind(w, 'radius');
            indxS = indxS(end) + 7;
            indxE = strfind(w(indxS:end), ' mm');
            indxE = indxE(1) - 2;
            SRad = w(indxS:indxS+indxE);
            
            % We don't extract the centre of gravity from here, since it needs
            % to be input in voxels... Instead get it from betted image
            Y = spm_read_vols(spm_vol(bSimg));
            Y = Y > 0;
            indY = find(Y);
            [subY_x subY_y subY_z] = ind2sub(size(Y), indY);
            COG = [mean(subY_x), mean(subY_y), mean(subY_z)];
            
            fprintf('\t...final c-o-g (vox): %0.4f %0.4f %0.4f  and radius (mm): %s\n', ...
                COG(1), COG(2), COG(3), SRad)
        end
        
        %% FIND OUTPUT
        % Get the mask images
        D = dir(fullfile(Spth, 'bet*mask*'));
        outMask = '';
        for d = 1:length(D)
            outMask = strvcat(outMask, fullfile(Spth, D(d).name));
        end
        
        % Get also the meshes
        D = dir(fullfile(Spth, 'bet*mesh*'));
        outMesh = '';
        for d = 1:length(D)
            outMesh = strvcat(outMesh, fullfile(Spth, D(d).name));
        end
        
        %% Make the BET BRAIN MASK
        % As made by BET [and slightly different from inskull_mask]
        V = spm_vol(bSimg);
        M = spm_read_vols(V);
        
        % Mask out non-brain
        M = M > 0;
        
        % Smooth mask if need be
        if aap.tasklist.currenttask.settings.smooth > 0
            fprintf('Smoothing the brain mask with %d voxel kernel \n', aap.tasklist.currenttask.settings.smooth)
            M = smooth3(M, 'box', repmat(aap.tasklist.currenttask.settings.smooth, [1 3]));
            M = M > 0;
        end
        
        % Then write out actual BET mask
        V.fname = fullfile(Spth, ['bet_' Sfn '_brain_mask' Sext]);
        spm_write_vol(V,M);
        % And add the mask to the list...
        outMask = strvcat(V.fname, outMask);
        
        %% MASK the brain(s)
        fprintf('Masking the brain(s) with Brain Mask \n')
        outStruct = '';
        for t = 1:length(aap.tasklist.currenttask.settings.structural)
            % Mask structural
            V = spm_vol(deblank(Simg(t,:)));
            Y = spm_read_vols(V);
            % Mask brain
            Y = Y.*M;
            % Write brain
            [pth nme ext]=fileparts(deblank(Simg(t,:)));
            V.fname = fullfile(pth,['bet_' nme ext]);
            spm_write_vol(V,Y);
            % Add to output...
            outStruct = strvcat(outStruct, V.fname);
        end
        
        %% DESCRIBE OUTPUTS!
        if aap.tasklist.currenttask.settings.maskBrain
            aap=aas_desc_outputs(aap,subj,inputstream{:},outStruct);
        else
            aap=aas_desc_outputs(aap,subj,inputstream{:},Simg);
        end
        maskStream = outputstream(~cellfun('isempty', strfind(outputstream,'BETmask')));
        aap=aas_desc_outputs(aap,subj, maskStream{:}, outMask);
        if aap.tasklist.currenttask.settings.masks
            meshStream = outputstream(~cellfun('isempty', strfind(outputstream,'BETmesh')));
            aap=aas_desc_outputs(aap,subj, meshStream{:}, outMesh);
        end
        
        %% DIAGNOSTIC IMAGE
        mriname = aas_prepare_diagnostic(aap,subj);
        
        %% Draw structural image...
        spm_check_registration(Simg)
        
        % This will only work for 1-7 masks
        OVERcolours = aas_colours;
        
        indx = 0;
        
        % Colour the brain extracted bit pink
        spm_orthviews('addcolouredimage',1,deblank(outStruct(1,:)), [0.9 0.4 0.4])
        % Add mesh outlines, to see if BET has worked properly!
        if aap.tasklist.currenttask.settings.masks
            for r = 1:size(outMesh,1)
                if strfind(outMesh(r,:), '.nii')
                    indx = indx + 1;
                    spm_orthviews('addcolouredimage',1,outMesh(r,:), OVERcolours{indx})
                end
            end
        else
            % Display outline of mask...
            copyfile(outMask(1,:), fullfile(Spth, 'betOutline.nii'));
            mask2outline(fullfile(Spth, 'betOutline.nii'));
            spm_orthviews('addcolouredimage', 1, fullfile(Spth, 'betOutline.nii'), OVERcolours{1})
        end
        
        spm_orthviews('reposition', [0 0 0])
        
        print('-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' mriname '.jpeg']));
        
        
        %% Diagnostic VIDEO of masks
        if aap.tasklist.currenttask.settings.diagnostic
            
            Ydims = {'X', 'Y', 'Z'};
            for d = 1:length(Ydims)
                aas_image_avi( Simg, ...
                    fullfile(pth, ['bet_' nme '_brain_mask' ext]), ...
                    fullfile(aap.acq_details.root, 'diagnostics', [mfilename '__' mriname '_' Ydims{d} '.avi']), ...
                    d, ... % Axis
                    [800 600], ...
                    2); % Rotations
            end
            try close(2); catch; end
        end
        
        % Clean up...
        if ~aap.tasklist.currenttask.settings.maskBrain
            delete(outStruct);
        end
end
