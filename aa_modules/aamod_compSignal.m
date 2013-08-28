% AA module wrapping...
% WARNING: This version requires you first to run the aamod_mask_fromstruct
% module present in the MVPaa toolbox and the aamod_fsl_BET module
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
    case 'domain'
        resp='session'; 
        
    case 'description'
        resp='Get signal from the CSF, WM, GM and OOB compartments';
        
    case 'summary'
        subjpath=aas_getsubjpath(subj);
        resp=sprintf('Get signal from the CSF, WM, GM and OOB compartments %s\n',subjpath);
        
    case 'report'
        
    case 'doit'
        
        % @@@ NOT YET IMPLEMENTED @@@ :
        % 1 - FOR OOB, WE WANT TOP CORNERS RELATIVE TO HEAD, TO AVOID GHOSTING
        
        % Let us use the native space...
        SMimg = aas_getfiles_bystream(aap,subj,'segmasksStrict');
        EPIimg = aas_getfiles_bystream(aap,subj,sess,'epi');
        BETimg = aas_getfiles_bystream(aap,subj,'epiBETmask');
        
        % Sanity checks
        [junk,fn] = fileparts(EPIimg(1,:));
        indx = strfind(fn, aap.directory_conventions.rawdataafterconversionprefix);
        indx = indx(1);
        if strfind(fn(1:indx-1), 'w')
            aas_log(aap, true, ['You should use unnormalised (i.e. native) images to do this analysis.' ...
                '\n\tThis should be run after the aamod_norm_noss and before aamod_norm_write'])
        end
        
        % Now, let's see which order the masks appear in...
        MOstr = aap.tasklist.currenttask.settings.maskOrder;
        MOlist = {};
        while ~isempty(MOstr)
            [tmp, MOstr] = strtok(MOstr,',');
            MOlist = [MOlist tmp];
        end
        
        % Load the segmented masks!
        mGM = []; mWM = []; mCSF = [];
        for m = 1:size(SMimg,1)
            % The ones *not* including string 'rwc' are the native masks
            [junk,fn] = fileparts(SMimg(m,:));
            if isempty(strfind(fn, 'rwc'))
                indx = strfind(fn, 'rc');
                eval(['m' MOlist{str2num(fn(indx + 2))} ' = spm_read_vols(spm_vol(SMimg(m,:)));'])
            end
        end
        
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
        
        % Record the number of voxels in each compartment
        nG = sum(mGM(:)>0);
        nW = sum(mWM(:)>0);
        nC = sum(mCSF(:)>0);
        nO = sum(mOOH(:)>0);
        
        fprintf('\nRemoving White Matter voxels near Gray Matter')
        mWM = rmNearVox(mWM, mGM, aap.tasklist.currenttask.settings.W2Gdist);
        
        % MASKS ALREADY STRICT ENOUGH TYPICALLY (TAKES OUT TOO MUCH CSF)
        %fprintf('\nRemoving CerebroSpinalFluid voxels near Gray Matter')
        %mCSF = rmNearVox(mCSF, mGM, aap.tasklist.currenttask.settings.C2Gdist);
        
        fprintf('\nRemoving CerebroSpinalFluid voxels near Skull')
        mCSF = rmNearVox(mCSF, mSkull, aap.tasklist.currenttask.settings.C2Sdist);
        
        fprintf('\nRemoving CerebroSpinalFluid voxels near OOH')
        mCSF = rmNearVox(mCSF, mOOH, aap.tasklist.currenttask.settings.C2Odist);
       
        %% Print the number of voxels in each compartment
        fprintf('\nGrey Matter mask comprises %d (%d) voxels', sum(mGM(:)>0), nG)
        fprintf('\nWhite Matter mask comprises %d (%d) voxels', sum(mWM(:)>0), nW)
        fprintf('\nCereberoSpinal Fluid mask comprises %d (%d) voxels', sum(mCSF(:)>0), nC)
        fprintf('\nOut of Head mask comprises %d (%d) voxels', sum(mOOH(:)>0), nO)
        
        compTC = zeros(size(EPIimg,1), 4);
        for e = 1:size(EPIimg,1)
            Y = spm_read_vols(spm_vol(EPIimg(e,:)));
            % Now average the data from each compartment
            compTC(e,1) = mean(Y(mGM>0));
            compTC(e,2) = mean(Y(mWM>0));
            compTC(e,3) = mean(Y(mCSF>0));
            compTC(e,4) = mean(Y(mOOH>0));
        end
        
        Rnames = {'GM', 'WM', 'CSF', 'OOH'};
        % Show an image of correlated timecourses...
        corrTCs(compTC, Rnames);

        %% DIAGNOSTIC IMAGE
        % Save graphical output to common diagnostics directory
        if ~exist(fullfile(aap.acq_details.root, 'diagnostics'), 'dir')
            mkdir(fullfile(aap.acq_details.root, 'diagnostics'))
        end
        mriname = strtok(aap.acq_details.subjects(subj).mriname, '/');
        set(gcf,'PaperPositionMode','auto')
        print('-djpeg','-r75',fullfile(aap.acq_details.root, 'diagnostics', ...
                [mfilename '__' mriname '.jpeg']));
            
            %% Diagnostic VIDEO of masks
        if aap.tasklist.currenttask.settings.diagnostic && ...
                sess == aap.acq_details.selected_sessions(end)
            
            movieFilename = fullfile(aap.acq_details.root, 'diagnostics', ...
                [mfilename '__' mriname '.avi']);
            % Create movie file by defining aviObject
            try delete(movieFilename); catch; end
            aviObject = avifile(movieFilename,'compression','none');
            
            mA = mGM + 2*mOOH + 3*mWM + 4*mCSF + 5*mSkull;
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
                aviObject = addframe(aviObject,getframe(2,windowSize));
            end

            aviObject = close(aviObject);
            try close(2); catch; end
        end
        
        %% DESCRIBE OUTPUTS!
        
        EPIdir = fileparts(EPIimg(1,:));
        save(fullfile(EPIdir, 'compSignal.mat'), 'compTC')
        aap=aas_desc_outputs(aap,subj,sess,'compSignal',fullfile(EPIdir, 'compSignal.mat'));
end
