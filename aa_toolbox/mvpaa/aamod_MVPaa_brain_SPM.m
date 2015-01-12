% AA module - Searchlight 1st level ==> THIS IS A FIRST LEVEL TASK
% You will need to make a local copy of this module into the same directory
%  as your user script, or somewhere else on your matlab path, and then modify it to
%  reflect the particular experiment you conducted
%
% Based on aa by Rhodri Cusack MRC CBU Mar 2006-Aug 2007
% Modified by Alejandro Vicente-Grabovetsky Dec-2008

function [aap,resp] = aamod_MVPaa_brain_SPM(aap,task,p)

resp='';

switch task
    case 'doit'
        
        fprintf('Working with data from participant %s. \n',aap.acq_details.subjects(p).mriname)
        
        Stats = []; EP = [];
        load(aas_getfiles_bystream(aap,p,'MVPaa'));
        
        % get sn mat file from normalisation
        normMAT = aas_getfiles_bystream(aap,p,'normalisation_seg_sn');
        
        % Load SPM used for this analysis...
        load(aas_getfiles_bystream(aap, p, 'firstlevel_spm'));
        
        % Example BET mask
        Mimg = aas_getfiles_bystream(aap, p, 'epiBETmask');
        % Brain mask!
        for a = 1:size(Mimg,1)
            if ~isempty(strfind(Mimg(a,:), 'brain_mask'))
                Mimg = deblank(Mimg(a,:));
                break
            end
        end           
        
        % FWHM in millimetres
        FWHMmm = aap.tasklist.currenttask.settings.FWHM;
        
        V = spm_vol(Mimg);
            
        brainSize = V.dim;
        
        %% GET MASK        
        mask = spm_read_vols(V);
        
        % Write out mask image...
        V.fname = fullfile(aas_getsubjpath(aap,p), 'mask.img');
        V.dt(1) = 2;
        spm_write_vol(V, mask);
        
        %% WRITE .img
        Stats = reshape(Stats, [brainSize(1), brainSize(2), brainSize(3), length(EP.contrasts), length(EP.tests)]);
        
        Flist = V.fname;
        V.dt(1) = 16; % Save in a format that accepts NaNs and negative values...
        
        fprintf('Saving images... \n')
        for c = 1:length(EP.contrasts)
            % Mean, median or beta
            V.fname = fullfile(aas_getsubjpath(aap,p), sprintf('con_%04d.img', c));
            Flist = strvcat(Flist, V.fname);
            spm_write_vol(V, squeeze(Stats(:,:,:,c,1)));
            
            % T-value
            V.fname = fullfile(aas_getsubjpath(aap,p), sprintf('spmT_%04d.img', c));
            Flist = strvcat(Flist, V.fname);
            spm_write_vol(V, squeeze(Stats(:,:,:,c,2)));
        end        
        
        %% NORMALISE
        if aap.tasklist.currenttask.settings.normalise == 1
            
            fprintf('Normalising images... \n')
            
            normPars = aap.spm.defaults.normalise.write;
            normPars.prefix = ''; % We want to keep no prefix...
            
            % This automatically reslices images to warped size
            spm_write_sn(Flist, normMAT, normPars);
        end
                
        %% SMOOTH IMAGES
        if FWHMmm > 0
            
            fprintf('Smoothing images... \n')
            for f = 2:size(Flist,1);
                Q = Flist(f,:);
                U = Flist(f,:); % No prefixes!
                spm_smooth(Q,U,FWHMmm);
            end
        end
        
        %% MASK SMOOTHED IMAGES!
        % Included mask to mask out untested data
        fprintf('NaNing untested voxels... \n')
        
        mask = spm_read_vols(spm_vol(Flist(1,:)));
        mask = mask > 0;
        
        for f = 2:size(Flist,1)
            V = spm_vol(Flist(f,:));
            Y = spm_read_vols(V);
            if strfind(Flist(f,:), 'spmT')
                % Zero mask in statistics...
                Y(~mask) = 0;
            elseif strfind(Flist(f,:), 'con')
                % NaN mask in statistics...
                Y(~mask) = NaN;
            end
            spm_write_vol(V, Y);
        end
                
        %% Modify SPM!
        % Clear SPM.xCon
        SPM.xCon = [];
        
        % Set correct path
        SPM.swd = aas_getsubjpath(aap,p);
        
        % Set world coordinates for visualisation...
        % ...which should already be found in the images...
        SPM.xVol.M = V.mat;
        SPM.xVol.iM = inv(SPM.xVol.M);
        
        % Size of the volume
        SPM.xVol.DIM = brainSize';
        
        % Smoothness of the volume...
        % ...Get the number of mm per voxel...
        mmVox = vox2mm(V);
        % ...then get the FWHM
        if FWHMmm < min(mmVox./2) % To avoid errors...
            FWHMmm = min(mmVox./2);
        end
        SPM.xVol.FWHM = [FWHMmm FWHMmm FWHMmm];
        SPM.xVol.FWHM = SPM.xVol.FWHM ./ mmVox;
        
        % Spm_resels_vol function
        % NOTE: This is probably not valid for FWE still, since the
        % searchlight procedure means each voxels is already "smoothed" to
        % some extent...
        SPM.xVol.R = spm_resels_vol( ...
            spm_vol(fullfile(aas_getsubjpath(aap,p), 'con_0001.img')), ...
            SPM.xVol.FWHM)';
        
        % Included voxels
        [X Y Z] = ind2sub(SPM.xVol.DIM',find(mask));
        SPM.xVol.XYZ = [X';Y';Z'];
        
        % Length of voxels in analysis
        SPM.xVol.S = length(X);
        
        % Filehandle of resels per voxel image (i.e. none!)
        SPM.xVol.VRpv = [];
                
        for c = 1:length(EP.contrasts)
            % SPM.xCon (.name)
            SPM.xCon(c).name = EP.contrasts(c).name;
            SPM.xCon(c).STAT = 'T';
            SPM.xCon(c).c = ones(size(SPM.xX.X,2),1);
            SPM.xCon(c).eidf = 1;
            SPM.xCon(c).Vcon = spm_vol(fullfile(aas_getsubjpath(aap,p), sprintf('con_%04d.img', c)));
            SPM.xCon(c).Vspm = spm_vol(fullfile(aas_getsubjpath(aap,p), sprintf('spmT_%04d.img', c)));
        end
        
        % Save SPM
        save(fullfile(aas_getsubjpath(aap,p), 'SPM.mat'), 'SPM');
        
        %% DESCRIBE OUTPUTS
        aap=aas_desc_outputs(aap,p,'firstlevel_spm', fullfile(aas_getsubjpath(aap,p), 'SPM.mat'));
        
        % Remove spmT images
        for f = size(Flist,1):-1:1
            if ~isempty(strfind(Flist(f,:), 'spmT_')) || ~isempty(strfind(Flist(f,:), 'mask'))
                Flist(f,:) = [];
            end
        end
        % Add headers to list of files...
        for f = 1:size(Flist,1)
            [Froot, Ffn, Fext] = fileparts(Flist(f,:));
            Flist = strvcat(Flist, fullfile(Froot, [Ffn '.hdr']));
        end
        aap=aas_desc_outputs(aap,p,'firstlevel_cons', Flist);
end