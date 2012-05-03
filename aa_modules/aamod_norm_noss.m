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
    case 'report'
        
    case 'doit'
        
        defs =aap.spm.defaults.normalise;
        defs.estimate.weight = '';
        
        % Template image; here template image not skull stripped
        % [AVG] Changed to allow specification of any T1 template, does not
        % need to be in the SPM folder any more...
        temp_imgs = aap.directory_conventions.T1template;
        if ~exist(temp_imgs, 'file')
            aas_log(aap, true, sprintf('Couldn''t find template T1 image %s.', temp_imgs));
        end
        
        clear imgs;
        
        %% Get structural
        % [AVG] Modified the way we get the structural, to be more aa4-like
        Simg = aas_getfiles_bystream(aap,subj,'structural');
        % Cheap and cheerful way of ensuring only one file is considered!
        if size(Simg,1) > 1
            Simg = deblank(Simg(1,:));
            aas_log(aap,0,sprintf('Found more than one structural so using:\n%s',Simg));
        end
        % Get structural directory for this subject
        [Spth, Sfn, Sext] = fileparts(Simg);
        
        %% Set up normalisation, etc.
        
        % Set the mask for subject to empty by default
        objMask = ''; % object mask
        
        % Because we are going reslice later (with undistort_reslice)
        % We don't reslice anything except the image to be normalized
        
        % call the SPM segment or normalize function to do the work
        if (aap.tasklist.currenttask.settings.usesegmentnotnormalise)
            % 2 stage process, as proposed by RH, to increased robustness [djm 13/03/06]
            
            %%%%%%%% 1st pass:
            fprintf('Running first pass of norm_noss (get bias corrected structural)\n')
            estopts.regtype='';    % turn off affine:
            out = spm_preproc(Simg, estopts);
            [sn,isn]   = spm_prep2sn(out);
            
            % only write out attenuation corrected image
            writeopts.biascor = 1;
            writeopts.GM  = [0 0 0];
            writeopts.WM  = [0 0 0];
            writeopts.CSF = [0 0 0];
            writeopts.cleanup = [0];
            spm_preproc_write(sn, writeopts);
            
            %%%%%%%% 2nd pass using attenuation corrected image
            fprintf('Running second pass of norm_noss (estimate normalisation)\n')
            % mstruc should be the attenuation corrected image
            % look for m prefixed filename
            mSimg = fullfile(Spth,['m' Sfn Sext]);
            if size(mSimg,1)>1
                aas_log(aap,0,sprintf('Found more than one attenuated structural so using first:\n%s',Simg));
            end
            
            estopts.regtype='mni';    % turn on affine again
            
            % Load header of image to be normalized
            V=spm_vol(mSimg);
            
            % Now adjust parameters according to starting offset parameters
            % as requested in task settings
            StartingParameters=[0 0 0   0 0 0   1 1 1   0 0 0];
            ParameterFields={'x','y','z', 'pitch','roll','yaw', 'xscale','yscale','zscale', 'xaffign','yaffign','zaffign'};
            fnames=fieldnames(aap.tasklist.currenttask.settings.affinestartingestimate);
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
            
            %[AVG] Save original V.mat parameters
            oldMAT = V.mat;
            
            % Adjust starting orientation of object image as requested
            StartingAffine=spm_matrix(StartingParameters);
            V.mat=StartingAffine*V.mat;
            
            % Run normalization
            out = spm_preproc(V,estopts);
            
            % Adjust output Affine to reflect fiddling of starting
            % orientation of object image
            out.Affine=out.Affine*StartingAffine;
            
            % [AVG] Instead we set the out.image parameters to our original
            % structural image!
            out.image.mat = oldMAT;
            
            [sn,isn]   = spm_prep2sn(out);
            
            % [AVG] DEBUG:
            % We could print out the Affines in an orderly way, so that
            % the experimenter can see what is being printed out...
            %{
            fprintf('\nInitial structural Affine\n')
            disp(oldMAT)
            fprintf('Initial Affine transform\n')
            disp(StartingAffine)
            fprintf('Out image Affine\n')
            disp(out.Affine)
            fprintf('Spatial Normalisation Affine\n')
            disp(sn.Affine)
            %}
            
            fprintf('Writing out the segmented images\n')
            % write out GM , WM, CSF native + unmod normalised
            writeopts.biascor = 1;
            writeopts.GM  = [0 1 1];
            writeopts.WM  = [0 1 1];
            writeopts.CSF = [0 1 1];
            writeopts.cleanup = [0];
            spm_preproc_write(sn,writeopts);
            
            SNmat = fullfile(Spth, [Sfn '_seg_sn.mat']);
            invSNmat = fullfile(Spth, [Sfn '_seg_inv_sn.mat']);
            savefields(SNmat,sn);
            savefields(invSNmat,isn);
            aap=aas_desc_outputs(aap,subj,'normalisation_seg_sn',SNmat);
            aap=aas_desc_outputs(aap,subj,'normalisation_seg_inv_sn',invSNmat);
            
            % [AVG] this output is completely different from .xml
            %{
            tiss={'grey','white','csf'};
            for tissind=1:3
                aap=aas_desc_outputs(aap,subj,sprintf('tissue_%s',tiss{tissind}),fullfile(Spth,sprintf('wc%d%s',tissind,mSimg)));
                aap=aas_desc_outputs(aap,subj,sprintf('unmod_tissue_%s',tiss{tissind}),fullfile(Spth,sprintf('c%d%s',tissind,mSimg)));
            end
            %}
            % [AVG] so instead, we group it all into segmentation stream
            outSeg = '';
            d = 0;
            while ~isnan(d)
                d = d+1;
                if exist(fullfile(Spth,sprintf('c%d%s',d,['m' Sfn Sext])), 'file')
                    outSeg = strvcat(outSeg, fullfile(Spth,sprintf('c%d%s',d,['m' Sfn Sext])));
                    outSeg = strvcat(outSeg, fullfile(Spth,sprintf('wc%d%s',d,['m' Sfn Sext])));
                else
                    d = NaN;
                end
            end
            aap=aas_desc_outputs(aap,subj,'segmentation',outSeg);
        else
            % Make the default normalization parameters file name
            % Turn off template weighting
            % SPM defaults
            SNmat = fullfile(Spth, [Sfn '_sn.mat']);
            invSNmat = fullfile(Spth, [Sfn '_seg_inv_sn.mat']);
            spm_normalise(temp_imgs, Simg, SNmat,...
                defs.estimate.weight, objMask, ...
                defs.estimate);
            aap=aas_desc_outputs(aap,subj,'normalisation_seg_sn',SNmat);
            aap=aas_desc_outputs(aap,subj,'normalisation_seg_inv_sn',invSNmat);
        end
        
        spm_write_sn(Simg,SNmat,defs.write);
        % [AVG] we need to add all the outputs, including warped structural
        % [AVG] It is probably best to save the 2ce bias-corrected image
        aap=aas_desc_outputs(aap,subj,'structural', strvcat( ...
            fullfile(Spth,['mm' Sfn Sext]), ...
            fullfile(Spth,['w' Sfn Sext])));
        
        %{
        % Now save graphical check
        figure(spm_figure('FindWin', 'Graphics'));
        % added graphical check for when segment is used [djm 20/01/06]
        if (aap.tasklist.currenttask.settings.usesegmentnotnormalise)
            myvols=spm_vol(char(aap.directory_conventions.T1template, ... % template T1
                Simg, ... % native T1
                fullfile(Spth,strcat('w',Simg)), ... % normalised T1
                fullfile(Spth,strcat('c1m',Simg)))); % native grey matter segmentation
            spm_check_registration(myvols)
            
            ann1=annotation('textbox',[.05 .96 .9 .03],'HorizontalAlignment','center','Color','r','String',strcat('Subject:...',Simg,',  processed on:...',date));
            ann2=annotation('textbox',[.1 .891 .3 .025],'HorizontalAlignment','center','Color','r','String','T1 template');
            ann3=annotation('textbox',[.6 .891 .3 .025],'HorizontalAlignment','center','Color','r','String','Native T1');
            ann4=annotation('textbox',[.1 .413 .3 .025],'HorizontalAlignment','center','Color','r','String','Normalised T1');
            ann5=annotation('textbox',[.6 .413 .3 .025],'HorizontalAlignment','center','Color','r','String','Native segmented grey matter');
            print('-djpeg',fullfile(subj_dir,'diagnostic_aamod_norm_noss'));
        end
        print('-djpeg',fullfile(subj_dir,'diagnostic_aamod_norm_noss'));
        if (aap.tasklist.currenttask.settings.usesegmentnotnormalise)
            delete(ann1); delete(ann2); delete(ann3);delete(ann4);delete(ann5);
        end
        %}
        
        % Diagnostic image?
        % Save graphical output to common diagnostics directory
        if ~exist(fullfile(aap.acq_details.root, 'diagnostics'), 'dir')
            mkdir(fullfile(aap.acq_details.root, 'diagnostics'))
        end
        mriname = strtok(aap.acq_details.subjects(subj).mriname, '/');
         try
            % This will only work for 1-7 segmentations
            OVERcolours = {[1 0 0], [0 1 0], [0 0 1], ...
                [1 1 0], [1 0 1], [0 1 1], [1 1 1]};
            
            %% Draw native template
            spm_check_registration(fullfile(Spth,['mm' Sfn Sext]))
            % Add normalised segmentations...
            for r = 1:(size(outSeg,1)/2)
                spm_orthviews('addcolouredimage',1,fullfile(Spth,sprintf('c%d%s',r, ['m' Sfn Sext])), OVERcolours{r})
            end
            %% Diagnostic VIDEO of segmentations
            aas_checkreg_avi(aap, subj, 2)
            
            spm_orthviews('reposition', [0 0 0])
            
            figure(spm_figure('FindWin', 'Graphics'));
            print('-djpeg','-r75',fullfile(aap.acq_details.root, 'diagnostics', ...
                [mfilename '__' mriname '_N.jpeg']));
            
            %% Draw warped template
            spm_check_registration(aap.directory_conventions.T1template)
            % Add normalised segmentations...
            for r = 1:(size(outSeg,1)/2)
                spm_orthviews('addcolouredimage',1,fullfile(Spth,sprintf('wc%d%s',r,['m' Sfn Sext])), OVERcolours{r})
            end
            spm_orthviews('reposition', [0 0 0])
            
            figure(spm_figure('FindWin', 'Graphics'));
            print('-djpeg','-r75',fullfile(aap.acq_details.root, 'diagnostics', ...
                [mfilename '__' mriname '_W.jpeg']));
        catch
            fprintf('\n\tFailed display diagnostic image - Displaying template & segmentation 1');
            try
                %% Draw native template
                spm_check_registration(char({fullfile(Spth,sprintf('c1%s',['m' Sfn Sext])); fullfile(Spth,['mm' Sfn Sext])}))
                
                figure(spm_figure('FindWin', 'Graphics'));
                print('-djpeg','-r75',fullfile(aap.acq_details.root, 'diagnostics', ...
                    [mfilename '__' mriname '_N.jpeg']));
                
                %% Draw warped template
                spm_check_registration(char({fullfile(Spth,sprintf('wc1%s',['m' Sfn Sext])); aap.directory_conventions.T1template}))
                
                figure(spm_figure('FindWin', 'Graphics'));
                set(gcf,'PaperPositionMode','auto')
                print('-djpeg','-r75',fullfile(aap.acq_details.root, 'diagnostics', ...
                    [mfilename '__' mriname ' W.jpeg']));
            catch
                fprintf('\n\tFailed display backup diagnostic image!');
            end
        end
        
        time_elapsed
        
    case 'checkrequirements'
        % Template image; here template image not skull stripped
        T1file = aap.directory_conventions.T1template;
        
        if (~exist(T1file,'file'))
            aas_log(aap,true,sprintf('T1 template file %s not found, check aap.directory_conventions.T1template\n',T1file));
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

return;
%------------------------------------------------------------------------