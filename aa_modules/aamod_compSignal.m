% AA module wrapping...
% WARNING: This version requires you first to run the
% aamod_mask_fromsegment and the aamod_fsl_BET module (optional)
%
%--------------------------------------------------------------------------
% COMP_SIGNAL creates regressors with mean signal intensity values for each
% segmented compartment [WhiteMatter (WM), CerebralSpinalFluid (CSF) and
% Out-of-Brain (OOB)] separately for each image. The resulting three
% regressors are saved in a mat file (structure: [WM CSF OOB]). These
% regressors can be used similarly to head motion regressors. Where the
% later can account for head motion effects, the compartment signal
% regressors can be used to account for global signal noise due to changes
% in the magnetic field over time (due to the movement of a conductive body
% - like an arm or hand - through the magnetic field) or other nuisance
% factors. Compartment signals are preferred over global signals as
% inclusion of the later might induce fake BOLD deactivations and a
% reduction in power. The compartment signals do not contain GreyMatter (so
% also no BOLD response) and therefore do not suffer from these ill
% effects.
%
% When using these regressors, you could cite my HBM abstract:
%   Verhagen L, Grol MJ, Dijkerman HC, Toni I. (2006) Studying visually-
%       guided reach to grasp movements in an MR-environment. Human Brain
%       Mapping.
%
% or cite my research paper which describes the methods superficially:
%   Verhagen L, Dijkerman HC, Grol MJ, Toni I (2008) Perceptuo-motor
%       interactions during prehension movements. J Neurosci 28(18):4726-4735
%
% or wait for the upcoming methods paper (a little patience required):
%   Verhagen L, Grol MJ, Dijkerman HC, Toni I. Studying visually-guided
%       reach to grasp movements in an MR-environment. Manuscript in
%       preparation.
%
% modified from comp_signal by Lennart Verhagen, 2005-01-25
% L.Verhagen@fcdonders.ru.nl
%--------------------------------------------------------------------------

function [aap,resp]=aamod_compSignal(aap,task,subj,sess)

resp='';

switch task
    
    case 'report'
        
    case 'doit'
        
        % @@@ NOT YET IMPLEMENTED @@@ :
        % 1 - FOR OOB, WE WANT TOP CORNERS RELATIVE TO HEAD, TO AVOID GHOSTING
        
        % Let us use the native space...
		inStreams = aas_getstreams(aap,'input');
        EPIimg = aas_getfiles_bystream(aap,subj,sess,'epi');
		SMimg = aas_getfiles_bystream(aap,subj,inStreams{2});        
        hasBET = false;
        if aas_stream_has_contents(aap,subj,'epiBETmask')
            BETimg = aas_getfiles_bystream(aap,subj,'epiBETmask');
            hasBET = true;
        end
        
        % Now, let's see which order the masks appear in...
        MOlist = textscan(aap.tasklist.currenttask.settings.maskOrder,'%s','delimiter',':');
        MOlist = MOlist{1}';
        
        % Load the segmented masks!
        mGM = []; mWM = []; mCSF = [];
        for m = 1:size(SMimg,1)
            [junk,fn] = fileparts(SMimg(m,:));
            indx = find(fn=='c',1,'first');
            feval(@()assignin('caller',['m' MOlist{str2num(fn(indx + 1))}],spm_read_vols(spm_vol(SMimg(m,:)))));
        end
        
        % Record the number of voxels in each compartment
        nG = sum(mGM(:)>0);
        nW = sum(mWM(:)>0);
        nC = sum(mCSF(:)>0);
        
        fprintf('\nRemoving White Matter voxels near Gray Matter\n')
        mWM = rmNearVox(mWM, mGM, aap.tasklist.currenttask.settings.W2Gdist);
        
        % MASKS ALREADY STRICT ENOUGH TYPICALLY (TAKES OUT TOO MUCH CSF)
        %fprintf('\nRemoving CerebroSpinalFluid voxels near Gray Matter')
        %mCSF = rmNearVox(mCSF, mGM, aap.tasklist.currenttask.settings.C2Gdist);
        
        if hasBET
            % Try to load the BET masks
            for m = 1:size(BETimg,1)
                [junk,fn] = fileparts(BETimg(m,:));
                if strfind(fn, 'outskin_mask')
                    mOOH = spm_read_vols(spm_vol(BETimg(m,:)));
                    mOOH = ~mOOH; % BET MASK IS INCLUSIVE HEAD...
                elseif strfind(fn, 'skull_mask')
                    mSkull = spm_read_vols(spm_vol(BETimg(m,:)));
                end
            end
            nO = sum(mOOH(:)>0);
            
            fprintf('Removing CerebroSpinalFluid voxels near Skull\n')
            mCSF = rmNearVox(mCSF, mSkull, aap.tasklist.currenttask.settings.C2Sdist);
            
            fprintf('Removing CerebroSpinalFluid voxels near OOH\n')
            mCSF = rmNearVox(mCSF, mOOH, aap.tasklist.currenttask.settings.C2Odist);
        end
        
        %% Print the number of voxels in each compartment
        fprintf('Grey Matter mask comprises %d (%d) voxels\n', sum(mGM(:)>0), nG)
        fprintf('White Matter mask comprises %d (%d) voxels\n', sum(mWM(:)>0), nW)
        fprintf('CereberoSpinal Fluid mask comprises %d (%d) voxels\n', sum(mCSF(:)>0), nC)
        if hasBET, fprintf('Out of Head mask comprises %d (%d) voxels\n', sum(mOOH(:)>0), nO); end
        
        if isfield(aap.options, 'NIFTI4D') && aap.options.NIFTI4D
            V = spm_vol(EPIimg);
        else
            for e = 1:size(EPIimg,1)
                V(e) = spm_vol(EPIimg(e,:));
            end
        end
        compTC = zeros(numel(V), 3);
        for e =  1:numel(V)
            Y = spm_read_vols(V(e));
            % Now average the data from each compartment
            compTC(e,1) = mean(Y(mGM>0));
            compTC(e,2) = mean(Y(mWM>0));
            compTC(e,3) = mean(Y(mCSF>0));
            if hasBET, compTC(e,4) = mean(Y(mOOH>0)); end
        end
        
        %% DESCRIBE OUTPUTS!
        sessdir = aas_getsesspath(aap,subj,sess);
        save(fullfile(sessdir, 'compSignal.mat'), 'compTC')
        aap=aas_desc_outputs(aap,subj,sess,'compSignal',fullfile(sessdir, 'compSignal.mat'));
        
        %% DIAGNOSTIC IMAGE
        Rnames = {'GM', 'WM', 'CSF'};
        if hasBET, Rnames{4} = 'OOH'; end
        
        % Show an image of correlated timecourses...
        corrTCs(compTC, Rnames);
        set(gcf,'PaperPositionMode','auto','Renderer','zbuffer');
        print('-djpeg','-r75',fullfile(sessdir, 'diagnostic_compSignal.jpg'));
        close(gcf);
        
        %% Diagnostic VIDEO of masks
        if strcmp(aap.options.wheretoprocess,'localsingle') && ...
                aap.options.diagnostic_videos && ...
                sess == aap.acq_details.selected_sessions(end)
            movieFilename = fullfile(sessdir, 'diagnostic_compSignal.avi');
            % Create movie file by defining aviObject
            try delete(movieFilename); catch; end
            if checkmatlabreq([7;11]) % From Matlab 7.11 use VideoWriter
                aviObject = VideoWriter(movieFilename);
                open(aviObject);
            else
                aviObject = avifile(movieFilename,'compression','none');
            end
            
            mA = mGM + 3*mWM + 4*mCSF;
            if hasBET
                mA = mA + 2*mOOH + 5*mSkull;
            end
            
            try close(2); catch; end
            figure(2)
            set(2, 'Position', [0 0 1000 800])
            windowSize = get(2,'Position');
            
            for x = 1:size(mA,1)
                h = subplot(1,1,1);
                imagesc(rot90(squeeze(mA(x,:,:))))
                axis equal off
                caxis([0 5])
                colorbar( ...
                    'Ytick', 0:5, ...
                    'Yticklabel', {'N/A', 'GM', 'OOB', 'WM', 'CSF', 'Skull'})
                zoomSubplot(h, 1.2)
                
                pause(0.01)
                % Capture frame and store in aviObject
                if checkmatlabreq([7;11]) % From Matlab 7.11 use VideoWriter
                    writeVideo(aviObject,getframe(2,windowSize));
                else
                    aviObject = addframe(aviObject,getframe(2,windowSize));
                end
            end
            if checkmatlabreq([7;11]) % From Matlab 7.11 use VideoWriter
                close(aviObject);
            else
                junk = close(aviObject);
            end
            try close(2); catch; end
        end
end
end