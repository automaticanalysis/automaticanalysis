% AA module
% Use the Anatomical Transformation Toolbox to normalise the structural to
% a template image

function [aap,resp]=aamod_ANTS_warp(aap,task,subj,sess)

resp='';

switch task
    case 'doit'
        subjname = aas_prepare_diagnostic(aap,subj);
        
        %% Get ANTs warps
        ANTSimg = aas_getfiles_bystream(aap,subj,'ANTs');
        % Get structural directory for this subject
        Spth = fileparts(ANTSimg(1,:));
        
        % Set the ANTS path
        setenv('ANTSPATH', aap.directory_conventions.ANTSdir)
        warpANTSpath = [fullfile(getenv('ANTSPATH'), 'bin', 'WarpImageMultiTransform') ' '];
        
        % Dimension number (always 3 for structural)
        Ndim = [num2str(3) ' '];
        
        if exist(fullfile(Spth,'antsAffine.txt'), 'file')
            affineTrans = [' ' fullfile(Spth, 'antsAffine.txt')]; % and affine, if this exists...
        else
            affineTrans = '';
        end
        
        %% Get images to warp!
        
        % The streams that we should normalise should be in the output...
        Ostreams = aap.tasklist.currenttask.outputstreams.stream;
        Istreams = aap.tasklist.currenttask.inputstreams.stream;
        
        % find out what stream we should use as a reference image
        Rstreams = Istreams(~strcmp('ANTs',Istreams));
        Rstreams = Rstreams(~strcmp(Ostreams,Rstreams));
        
        if length(Rstreams) > 1
            aas_log(aap,true', 'Too many reference images')
        end
        
        %% Reference image to which your warped images are resliced to!
        if isempty(Rstreams)
            % If there's no reference stream, then use template...
            refimg = aap.directory_conventions.T1template;
        else
            refimg = aas_getfiles_bystream(aap,Rstreams{:});
        end
        % Use first image of stream as reference...
        refimg = refimg(1,:);
        
        if ~exist(refimg, 'file')
            aas_log(aap, true, sprintf('Couldn''t find reference image %s.', refimg));
        end
        
        for streamind=1:length(Ostreams)
            
            % Images to warp
            if (exist('sess','var'))
                P = aas_getfiles_bystream(aap,subj,sess,Ostreams{streamind});
            else
                P = aas_getfiles_bystream(aap,subj,Ostreams{streamind});
            end
            
            wimgs={};
            
            for c=1:size(P,1)
                [pth fn ext]=fileparts(deblank(P(c,:)));
                wimgs = [wimgs,fullfile(pth,['w' fn ext])]; % [TNS] name of the output file ? why do we have to delete a non-existing file afterwards?
                
                % delete previous because otherwise nifti write routine doesn't
                % save disc space when you reslice to a coarser voxel
                [s w]=aas_shell(['rm ' fullfile(pth,['w' fn ext])],true); % quietly
                
                %% Use ANTS to normalise the stream!
                
                if ~aap.tasklist.currenttask.settings.inversetransform % [TNS] added inverse transform option
                    warpANTS_command = [ warpANTSpath Ndim ... % dimension number
                        P(c,:) ' ' fullfile(pth, ['w' fn ext]) ... % moving image & output
                        ' -R ' refimg ' '... % reference image is the template
                        fullfile(Spth, 'antsWarp.nii') ' '... % affine transform (if it exists)
                        affineTrans ]; % transform
                else
                    warpANTS_command = [ warpANTSpath Ndim ... % dimension number
                        P(c,:) ' ' fullfile(pth, ['w' fn ext]) ... % moving image & output
                        ' -R ' refimg ' '... % reference image now is the subjects' structural image
                        '-i ' affineTrans ' '... % inverse affine transform (if it exists)
                        fullfile(Spth, 'antsInverseWarp.nii')]; % transform
                end
                [s w] = aas_shell(warpANTS_command);
            end
            
            if isempty(aap.tasklist.currenttask.settings.checkreg)
                aap.tasklist.currenttask.settings.checkreg = 1:length(wimgs);
            end
            %% describe outputs
            
            if (exist('sess','var'))
                aap=aas_desc_outputs(aap,subj,sess,Ostreams{streamind},wimgs);
            else
                aap=aas_desc_outputs(aap,subj,Ostreams{streamind},wimgs);
            end
            
            warning off
            for c = aap.tasklist.currenttask.settings.checkreg
                %% Draw native template
                spm_check_registration(strvcat( ...
                    wimgs{c}, ...
                    refimg))
                
                %% Diagnostic IMAGE of segmentations
                
                % Save graphical output to common diagnostics directory
                if ~exist(fullfile(aap.acq_details.root, 'diagnostics'), 'dir')
                    mkdir(fullfile(aap.acq_details.root, 'diagnostics'))
                end
                
                spm_orthviews('reposition', [0 0 0])
                spm_ov_reorient('context_init', 1)
                
                try figure(spm_figure('FindWin', 'Graphics')); catch; figure(1); end;
                print('-depsc2', fullfile(aap.acq_details.root, 'diagnostics', ...
                    [mfilename '__' subjname '_' num2str(c) '.eps']));
            end
            warning on
        end
end
end