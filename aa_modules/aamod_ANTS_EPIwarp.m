% AA module
% Use the Anatomical Transformation Toolbox to normalise the structural to
% a template image

function [aap,resp]=aamod_ANTS_EPIwarp(aap,task,subj,sess)

resp='';

switch task
    case 'doit'
        
        warning off
        
        % Set the ANTS path
        setenv('ANTSPATH', aap.directory_conventions.ANTSdir)
        warpANTSpath = [fullfile(getenv('ANTSPATH'), 'bin', 'WarpImageMultiTransform') ' '];
        
        % Dimension number (always 3 for structural)
        Ndim = [num2str(3) ' '];
        
        %% Get ANTs warps       
        ANTSimg = aas_getfiles_bystream(aap,subj,'ANTs');
        % Get structural directory for this subject
        Spth = fileparts(ANTSimg(1,:));
        
        %% Template image
        % [AVG] Changed to allow specification of any T1 template, does not
        % need to be in the SPM folder any more...
        sTimg = aap.directory_conventions.T1template;
        if ~exist(sTimg, 'file')
            aas_log(aap, true, sprintf('Couldn''t find template T1 image %s.', Timg));
        end
        
        %% Get EPI images (or meanEPI, or...)
        
        streams=aap.tasklist.currenttask.inputstreams;
        
        % find out what streams we should normalise
        streams=streams.stream(~[strcmp('ANTs',streams.stream)]);
        
        % Is session specified in task header (used for meanepi, which only
        % occurs in session 1)
        if (isfield(aap.tasklist.currenttask.settings,'session'))
            sess=aap.tasklist.currenttask.settings.session;
        end
        
        for streamind=1:length(streams)
            
            % Images to warp
            if (exist('sess','var'))
                P = aas_getfiles_bystream(aap,subj,sess,streams{streamind});
            else
                P = aas_getfiles_bystream(aap,subj,streams{streamind});
            end;
            
            wimgs=[];
            
            for c=1:size(P,1)
                [pth fn ext]=fileparts(P(c,:));
                wimgs = strvcat(wimgs,fullfile(pth,['w' fn ext]));
                
                % delete previous because otherwise nifti write routine doesn't
                % save disc space when you reslice to a coarser voxel
                [s w]=aas_shell(['rm ' fullfile(pth,['w' fn ext])],true); % quietly
                
                %% Use ANTS to normalise the stream!
                warpANTS_command = [ warpANTSpath Ndim ... % dimension number
                    P(c,:) ' ' fullfile(pth, ['w' fn ext]) ... % moving image & output
                    ' -R ' sTimg ' '... % reference image
                    fullfile(Spth, 'antsWarp.nii')]; % transform
                if exist(fullfile(Spth,'antsAffine.txt'), 'file')
                    warpANTS_command = [warpANTS_command ' ' fullfile(Spth, 'antsAffine.txt')]; % and affine, if this exists...
                end
                
                [s w] = aas_shell(warpANTS_command);
            end
            
            %% describe outputs
            
            if (exist('sess','var'))
                aap=aas_desc_outputs(aap,subj,sess,streams{streamind},wimgs);
            else
                aap=aas_desc_outputs(aap,subj,streams{streamind},wimgs);
            end
        end
end
end