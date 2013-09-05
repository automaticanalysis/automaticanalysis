% AA module - normalisation using normalise or two pass procedure with segment
% [aap,resp]=aamod_norm_noss(aap,task,subj)
% Depending on aap.tasksettings.aamod_norm_noss.usesegmentnotnormalise
% If 0 - classic Normalisation
% If 1 - use Segment with two pass procedure: a first pass to correct
%  for intensity variation across our structural images (probably due to
%  inhomogeneous SNR across space of 12 channel coil); and a second pass
%  to then do the segmentation
% _noss version does not use skull stripping
% subj=subject num
% Rhodri Cusack & Daniel Mitchell MRC CBU 2006
% based on originals by Rik Henson, Matthew Brett

function [aap,resp]=aamod_norm_noss(aap,task,subj)
resp='';

switch task
    case 'report' % [TA]
        if ~exist(fullfile(aas_getsubjpath(aap,subj),['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '.jpg']),'file')
            diag(aap,subj);
        end
        fdiag = dir(fullfile(aas_getsubjpath(aap,subj),'diagnostic_*.jpg'));
        for d = 1:numel(fdiag)
            aap = aas_report_add(aap,subj,'<table><tr><td>');
            aap=aas_report_addimage(aap,subj,fullfile(aas_getsubjpath(aap,subj),fdiag(d).name));
            aap = aas_report_add(aap,subj,'</td></tr></table>');
        end
    case 'doit'
        
        defs =aap.spm.defaults.normalise;
        defs.estimate.weight = '';
        
        % Template image
        temp_imgs = aap.directory_conventions.T1template;
        if (~exist(temp_imgs,'file')) % try in SPM
            if temp_imgs(1) ~= '/', temp_imgs = fullfile(fileparts(which('spm')),temp_imgs); end
        end
        if ~exist(temp_imgs, 'file')
            aas_log(aap, true, sprintf('Couldn''t find template T1 image %s.', temp_imgs));
        end
        
        clear imgs;
        
        %% Get input structural image
        streams=aap.tasklist.currenttask.inputstreams.stream{:};
        
        Simg = aas_getfiles_bystream(aap,subj,streams);
        % Which file is considered, as determined by the structural parameter!
        if size(Simg,1) > 1
            Simg = deblank(Simg(aap.tasklist.currenttask.settings.structural, :));
            fprintf('WARNING: Several structurals found, considering: \n')
            fprintf('\t%s\n', Simg(1,:))
        end
        % Get structural directory for this subject
        [Spth, Sfn, Sext] = fileparts(Simg);
        
        % We will be getting SN transformation matrices...
        SNmat = fullfile(Spth, [Sfn '_seg_sn.mat']);
        invSNmat = fullfile(Spth, [Sfn '_seg_inv_sn.mat']);
        
        %% Set up normalisation, etc.
        
        % Set the mask for subject to empty by default
        objMask = ''; % object mask
        
        % Because we are going reslice later (with undistort_reslice)
        % We don't reslice anything except the image to be normalized
        
        % call the SPM segment or normalize function to do the work
        if aap.tasklist.currenttask.settings.usesegmentnotnormalise
            if aap.tasklist.currenttask.settings.bias1stpass
                % 2 stage process, as proposed by RH, to increased robustness
                % [djm 13/03/06]
                %%%%%%%% 1st pass:
                fprintf('Running first pass of norm_noss (get bias corrected structural)\n')
                estopts.regtype='';    % turn off affine:
                if ~isfield(estopts,'fudge'), estopts.fudge = 5; end % compatiblity
                out = spm_preproc(Simg, estopts);
                [sn,isn]   = spm_prep2sn(out);
                
                % only write out attenuation corrected image
                writeopts.biascor = 1;
                writeopts.GM  = [0 0 0];
                writeopts.WM  = [0 0 0];
                writeopts.CSF = [0 0 0];
                writeopts.cleanup = [0];
                spm_preproc_write(sn, writeopts);
                
                mSimg = fullfile(Spth,['m' Sfn Sext]);
            else
                mSimg = Simg;
            end
            [mSpth, mSfn, mSext] = fileparts(mSimg);
            
            %%%%%%%% 2nd pass using attenuation corrected image
            fprintf('Running second pass of norm_noss (estimate normalisation)\n')
            % mstruc should be the attenuation corrected image
            % look for m prefixed filename
            
            estopts.regtype='mni';    % turn on affine again
            
            % Load header of image to be normalized
            V = spm_vol(mSimg);
            
            % Now adjust parameters according to starting offset parameters as requested in task settings
            StartingParameters=[0 0 0   0 0 0   1 1 1   0 0 0];
            ParameterFields={'x','y','z', 'pitch','roll','yaw', 'xscale','yscale','zscale', 'xaffign','yaffign','zaffign'};
            if ~isempty(aap.tasklist.currenttask.settings.affinestartingestimate)
                fnames=fieldnames(aap.tasklist.currenttask.settings.affinestartingestimate);
            else
                fnames = [];
            end
            for fieldind=1:length(fnames)
                % Which element in StartingParameters does this refer to?
                whichitem=find([strcmp(fnames{fieldind},ParameterFields)]);
                % Generate a helpful error if it isn't recognised
                if (isempty(whichitem))
                    err=sprintf('Unexpected field %s in header file aamod_norm_noss.xml - expected one of ',fnames{fieldind});
                    err=[err sprintf('%s\t',ParameterFields)];
                    aas_log(aap,true,err);
                end
                % Put this in its place
                StartingParameters(whichitem)=aap.tasklist.currenttask.settings.affinestartingestimate.(fnames{fieldind});
            end
            
            % Save original V.mat parameters
            oldMAT = V.mat;
            
            % Adjust starting orientation of object image as requested
            StartingAffine=spm_matrix(StartingParameters);
            V.mat=StartingAffine*V.mat;
            
            % Run normalization
            out = spm_preproc(V,estopts);
            
            % Adjust output Affine to reflect fiddling of starting orientation of object image
            out.Affine=out.Affine*StartingAffine;
            
            % Instead we set the out.image parameters to our original structural image!
            out.image.mat = oldMAT;
            
            [sn,isn] = spm_prep2sn(out);
            
            % [AVG] DEBUG:
            % Let us print out the Affines in an orderly way, so that
            % the experimenter can see what is being printed out...
            fprintf('\nInitial structural Affine\n')
            disp(oldMAT)
            fprintf('Initial Affine transform\n')
            disp(StartingAffine)
            fprintf('Out image Affine\n')
            disp(out.Affine)
            fprintf('Spatial Normalisation Affine\n')
            disp(sn.Affine)
            
            fprintf('Writing out the segmented images\n')
            % write out GM, WM, CSF native + unmod normalised
            writeopts.biascor = 1;
            writeopts.GM  = [0 1 1];
            writeopts.WM  = [0 1 1];
            writeopts.CSF = [0 1 1];
            writeopts.cleanup = [0];
            spm_preproc_write(sn,writeopts);
            
            savefields(SNmat,sn);
            savefields(invSNmat,isn);
            
            % We group it all into segmentation stream
            outSeg = '';
            d = 0;
            while ~isnan(d)
                d = d+1;
                if exist(fullfile(mSpth,sprintf('c%d%s',d,[mSfn mSext])), 'file')
                    outSeg = strvcat(outSeg, fullfile(mSpth,sprintf('c%d%s',d,[mSfn mSext])));
                    outSeg = strvcat(outSeg, fullfile(mSpth,sprintf('wc%d%s',d,[mSfn mSext])));
                else
                    d = NaN;
                end
            end
            aap=aas_desc_outputs(aap,subj,'segmentation',outSeg);
        else
            % Make the default normalization parameters file name
            % Turn off template weighting
            % SPM defaults
            
            spm_normalise(temp_imgs, Simg, SNmat,...
                defs.estimate.weight, objMask, ...
                defs.estimate);
        end
        
        if ~all(cellfun('isempty', strfind(aap.tasklist.currenttask.outputstreams.stream, 'normalisation')))
            
            aap=aas_desc_outputs(aap,subj,'normalisation_seg_sn',SNmat);
            % SPM2 normalization doesn't generate the inverse transformation
            try
                aap=aas_desc_outputs(aap,subj,'normalisation_seg_inv_sn',invSNmat);
            catch
            end
            
            spm_write_sn(Simg,SNmat,defs.write);
            
            % [AVG] we need to add all the outputs, including warped structural
            % [AVG] It is probably best to also save the last bias-corrected image
            % [CW] But we don't have a bias corrected image if we didn't use
            % segmentation.
            if (aap.tasklist.currenttask.settings.usesegmentnotnormalise)
                aap=aas_desc_outputs(aap,subj,streams, strvcat( ...
                    fullfile(mSpth,['m' mSfn mSext]), ...
                    fullfile(Spth,['w' Sfn Sext])));
            else
                aap=aas_desc_outputs(aap,subj,streams, fullfile(Spth,['w' Sfn Sext]));
            end
            
            ann1=annotation('textbox',[.05 .96 .9 .03],'HorizontalAlignment','center','Color','r','String',strcat('Subject:...',Simg,',  processed on:...',date));
            ann2=annotation('textbox',[.1 .891 .3 .025],'HorizontalAlignment','center','Color','r','String','T1 template');
            ann3=annotation('textbox',[.6 .891 .3 .025],'HorizontalAlignment','center','Color','r','String','Native T1');
            ann4=annotation('textbox',[.1 .413 .3 .025],'HorizontalAlignment','center','Color','r','String','Normalised T1');
            ann5=annotation('textbox',[.6 .413 .3 .025],'HorizontalAlignment','center','Color','r','String','Native segmented grey matter');
            print('-djpeg',fullfile(Spth,'diagnostic_aamod_norm_noss'));
        end
        print('-djpeg',fullfile(Spth,'diagnostic_aamod_norm_noss'));
        if (aap.tasklist.currenttask.settings.usesegmentnotnormalise)
            delete(ann1); delete(ann2); delete(ann3);delete(ann4);delete(ann5);
        end
        
        %% Save graphical output to common diagnostics directory
        diag(aap,subj);
        
    case 'checkrequirements'
        % Template image; here template image not skull stripped
        T1file = aap.directory_conventions.T1template;
        if (~exist(T1file,'file')) % try in SPM
            if T1file(1) ~= '/', T1file = fullfile(fileparts(which('spm')),T1file); end
        end
        if (~exist(T1file,'file'))
            aas_log(aap,true,sprintf('T1 template file %s not found, check aap.directory_conventions.T1template\n',T1file));
        end
end
end
%------------------------------------------------------------------------
function savefields(fnam,p)
if length(p)>1, error('Can''t save fields.'); end
fn = fieldnames(p);
if numel(fn)==0, return; end
for subj=1:length(fn),
    eval([fn{subj} '= p.' fn{subj} ';']);
end
if str2double(version('-release'))>=14,
    save(fnam,'-V6',fn{:});
else
    save(fnam,fn{:});
end
end
%------------------------------------------------------------------------
function diag(aap,subj) % [TA]
mriname = aas_prepare_diagnostic(aap,subj);

% SPM, AA
Simg = aas_getfiles_bystream(aap,subj,'structural');
[Spth, Sfn, Sext] = fileparts(Simg);
outSeg = aas_getfiles_bystream(aap,subj,'segmentation');
try
    OVERcolours = aas_colours(size(outSeg,1)/2);
    
    %% Draw native template
    spm_check_registration(fullfile(mSpth,['m' mSfn mSext]))
    % Add segmentations...
    for t = 1:(size(outSeg,1)/2)
        spm_orthviews('addcolouredimage',1,fullfile(mSpth,sprintf('c%d%s',t, [mSfn mSext])), OVERcolours{t})
    end
    
    spm_orthviews('reposition', [0 0 0])
    
    try figure(spm_figure('FindWin', 'Graphics')); catch; figure(1); end;
    print('-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
        [mfilename '__' mriname '_N.jpeg']));
    
    %% Draw warped template
    spm_check_registration(aap.directory_conventions.T1template)
    % Add normalised segmentations...
    for t = 1:(size(outSeg,1)/2)
        spm_orthviews('addcolouredimage',1,fullfile(mSpth,sprintf('wc%d%s',t,[mSfn mSext])), OVERcolours{t})
    end
    spm_orthviews('reposition', [0 0 0])
    
    try figure(spm_figure('FindWin', 'Graphics')); catch; figure(1); end;
    print('-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
        [mfilename '__' mriname '_W.jpeg']));
catch
    fprintf('\n\tFailed display diagnostic image - Displaying template & segmentation 1');
    try
        %% Draw native template
        spm_check_registration(char({fullfile(Spth,sprintf('c1%s',['m' Sfn Sext])); ...
            fullfile(Spth,['mm' Sfn Sext])}))
        
        try figure(spm_figure('FindWin', 'Graphics')); catch; figure(1); end;
        print('-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' mriname '_N.jpeg']));
        
        %% Draw warped template
        spm_check_registration(char({fullfile(Spth,sprintf('wc1%s',['m' Sfn Sext])); ...
            aap.directory_conventions.T1template}))
        
        try figure(spm_figure('FindWin', 'Graphics')); catch; figure(1); end;
        set(gcf,'PaperPositionMode','auto')
        print('-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' mriname '_W.jpeg']));
    catch
        fprintf('\n\tFailed display backup diagnostic image!');
    end
end

% Another diagnostic image, looking at how well the segmentation worked...
if aap.tasklist.currenttask.settings.usesegmentnotnormalise
    Pthresh = 0.95;
    
    ROIdata = roi2hist(fullfile(Spth,['m' Sfn Sext]), ...
        outSeg(1:2:end,:), Pthresh);
    
    [h, pv, ci, stats] = ttest2(ROIdata{2}, ROIdata{1});
    
    title(sprintf('GM vs WM... T-val: %0.2f (df = %d)', stats.tstat, stats.df))
    
    print('-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
        [mfilename '__' mriname '_Hist.jpeg']));
end

%% Diagnostic VIDEO
if aap.tasklist.currenttask.settings.diagnostic
    Ydims = {'X', 'Y', 'Z'};
    
    for d = 1:length(Ydims)
        if (aap.tasklist.currenttask.settings.usesegmentnotnormalise)
            aas_image_avi(fullfile(mSpth,['m' mSfn mSext]), ...
                outSeg([1:2:size(outSeg,1)],:), ...
                fullfile(aap.acq_details.root, 'diagnostics', [mfilename '__' mriname '_' Ydims{d} '.avi']), ...
                d, ... % Axis
                [800 600], ...
                2, ... % Rotations
                'none'); % No outline...
            try close(2); catch; end
        end
    end
end

% FSL
% Obtain FSL T1 template
tP = fullfile(getenv('FSLDIR'),'data','standard','MNI152_T1_2mm_brain.nii.gz');

% Obtain normalized GM segmentation
Spth=aas_getsubjpath(aap,subj);
segdir=fullfile(Spth,aap.directory_conventions.structdirname);
sP = dir( fullfile(segdir,['wc1' aap.acq_details.subjects(subj).structuralfn '*.nii']));
sP = fullfile(segdir,sP(1).name);

% Create FSL-like overview
iP = fullfile(Spth,['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name]);
aas_runfslcommand(aap,sprintf('slices %s %s -s 2 -a %s.gif',tP,sP,iP));
[img,map] = imread([iP '.gif']); s3 = size(img,1)/3;
img = horzcat(img(1:s3,:,:),img(s3+1:2*s3,:,:),img(s3*2+1:end,:,:));
imwrite(img,[iP '.jpg']); delete([iP '.gif']);
end
