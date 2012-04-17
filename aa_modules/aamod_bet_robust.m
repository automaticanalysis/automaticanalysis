% AA module
% Runs BET (FSL Brain Extration Toolbox) on structural
% ROBUST:
% 1) Coregisters the T1 template to the structural and then cuts out neck
% from structural
% 2) Runs the -r option on BET to obtain a robust segmentation
% 3) Outputs brain masks and meshes
% [For correct functionality, it is recommended you run this after
% realignment and before writing the normalised EPI image
% If you do it before estimating the normalisation, make sure you normalise
% to a scull-stripped template, if at all possible!]

function [aap,resp]=aamod_bet_robust(aap,task,p)

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
        Simg = aas_getfiles_bystream(aap,p,'structural');
        mEPIimg = aas_getfiles_bystream(aap,p,1,'meanepi');
        
        % Which file is considered, as determined by the structural parameter!
        if size(Simg,1) > 1
            Simg = deblank(Simg(aap.tasklist.currenttask.settings.structural, :));
            fprintf('\tWARNING: Several structurals found, considering: %s\n', Simg)
        end
        % With the mean EPI, we just use the first one (there really should be only one)
        if size(mEPIimg,1) > 1
            mEPIimg = deblank(mEPIimg(1,:));
            fprintf('\tWARNING: Several mean EPIs found, considering: %s\n', mEPIimg)
        end
        
        if aap.tasklist.currenttask.settings.premask
            Sdir = fileparts(Simg);
            
            %% 0) Check that the tamplates we need exist!
            % Get the template
            sTimg = fullfile(spm('dir'), 'templates', 'T1.nii');
            if ~exist(sTimg, 'file')
                aas_log(aap, true, sprintf('Couldn''t find template T1 image %s.', Timg));
            end
            
            %% 1) Structural to T1 template
            fprintf('Coregistering the T1 template to structural\n')
            
            global defaults %#ok<TLEV>
            flags = defaults.coreg;
            
            % Copy template to structural location
            copyfile(sTimg, Sdir);
            sTimg = fullfile(Sdir, 'T1.nii');
            
            % Coregister template to Structural
            x = spm_coreg(spm_vol(Simg), spm_vol(sTimg), flags.estimate);
            
            % Set the new space for the template
            MM = spm_get_space(sTimg);
            spm_get_space(sTimg, spm_matrix(x)\MM);
            
            %% 2) Then reslice the Template
            
            fprintf('Reslicing T1 template to structural\n')
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
            
            % Reslice
            spm_reslice(strvcat(Simg, sTimg), resFlags);
            
            %% 3) Mask the Structural image using the T1 template
            
            fprintf('Mask structural with resliced T1 template\n')
            
            M = spm_read_vols(spm_vol(fullfile(Sdir, 'rT1.nii')));
            M = M>0;
            
            % Mask structural
            V = spm_vol(Simg);
            Y = spm_read_vols(V);
            Y = Y.*M;
            spm_write_vol(V,Y);
            
            delete(fullfile(Sdir, 'T1.nii'));
            delete(fullfile(Sdir, 'rT1.nii'));
        end
        
        [pth nme ext]=fileparts(Simg);
        
        outStruct=fullfile(pth,['bet_' nme ext]);
        % Run BET [-R Using robust setting to avoid neck!]   
        fprintf('1st BET pass (recursive) to find optimal centre of gravity and radius\n')
        [junk, w]=aas_runfslcommand(aap, ...
            sprintf('bet %s %s -f %f -v -R',Simg,outStruct, ...
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
        
        if aap.tasklist.currenttask.settings.masks
            fprintf('2nd BET pass extracting brain masks \n')
            % Run BET [-A Now let's get the brain masks and meshes!!]
            [junk, w]=aas_runfslcommand(aap, ...
                sprintf('bet %s %s -f %f -c %0.4f %0.4f %0.4f -r %s -v -A',Simg,outStruct, ...
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
            Y = spm_read_vols(spm_vol(outStruct));
            Y = Y > 0;
            indY = find(Y);
            [subY_x subY_y subY_z] = ind2sub(size(Y), indY);
            COG = [mean(subY_x), mean(subY_y), mean(subY_z)];
            
            fprintf('\t...final c-o-g (vox): %0.4f %0.4f %0.4f  and radius (mm): %s\n', ...
                COG(1), COG(2), COG(3), SRad)
        end
        
        % BRAIN MASK (slightly different from inskull)
        
        V = spm_vol(fullfile(pth,['bet_' nme ext]));
        mY = spm_read_vols(V);
        
        % Mask out non-brain
        mY = mY > 0;
        
        % Then write out actual BET mask
        V.fname = fullfile(pth, ['bet_' nme '_brain_mask' ext]);
        spm_write_vol(V,mY);
        
        
        %% DIAGNOSTIC IMAGE
        % Save graphical output to common diagnostics directory
        if ~exist(fullfile(aap.acq_details.root, 'diagnostics'), 'dir')
            mkdir(fullfile(aap.acq_details.root, 'diagnostics'))
        end
        mriname = strtok(aap.acq_details.subjects(p).mriname, '/');
        try
            %% Draw structural image...
            spm_check_registration(Simg)
            
            % This will only work for 1-7 masks
            OVERcolours = {[1 0 0], [0 1 0], [0 0 1], ...
                [1 1 0], [1 0 1], [0 1 1], [1 1 1]};
            
            indx = 0;
            
            % Colour the brain extracted bit pink
            spm_orthviews('addcolouredimage',1,outStruct, [0.9 0.4 0.4])
            % Add mesh outlines, to see if BET has worked properly!
            if aap.tasklist.currenttask.settings.masks
                for r = 1:size(outMesh,1)
                    if strfind(outMesh(r,:), '.nii')
                        indx = indx + 1;
                        spm_orthviews('addcolouredimage',1,outMesh(r,:), OVERcolours{indx})
                    end
                end
            end
            %% Diagnostic VIDEO of masks
            aas_checkreg_avi(aap, p, 2)
            
            spm_orthviews('reposition', [0 0 0])
            
            figure(spm_figure('FindWin'));
            set(gcf,'PaperPositionMode','auto')
            print('-djpeg','-r75',fullfile(aap.acq_details.root, 'diagnostics', ...
                [mfilename '__' mriname '.jpeg']));
        catch
        end
        
        %% DESCRIBE OUTPUTS!
        
        % Get the mask images
        D = dir(fullfile(pth, 'bet*mask*'));
        outMask = '';
        for d = 1:length(D)
            outMask = strvcat(outMask, fullfile(pth, D(d).name));
        end
        
        % Get also the meshes
        D = dir(fullfile(pth, 'bet*mesh*'));
        outMesh = '';
        for d = 1:length(D)
            outMesh = strvcat(outMesh, fullfile(pth, D(d).name));
        end
        
        % Structural image after BETting
        aap=aas_desc_outputs(aap,p,'structural',outStruct);
        aap=aas_desc_outputs(aap,p,'BETmask',outMask);
        if aap.tasklist.currenttask.settings.masks
            aap=aas_desc_outputs(aap,p,'BETmesh',outMesh);
        end
        
        time_elapsed
end
