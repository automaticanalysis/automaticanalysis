% AA module
% Convert the fieldmap images (2 mag and 1 phase) into a Voxel Displacement
% Map (VDM) using the FieldMap toolbox of SPM

function [aap,resp]=aamod_fieldmap2VDM(aap,task,p)

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
        
        % Defaults specified in this path
        % You can set your own settings in your own copy of the XML or recipe!        
        pm_defs = ...
            [aap.tasklist.currenttask.settings.te1, ... % te1
            aap.tasklist.currenttask.settings.te2, ... % te2
            aap.tasklist.currenttask.settings.epifm, ... % epifm
            aap.tasklist.currenttask.settings.tert, ... % tert
            aap.tasklist.currenttask.settings.kdir, ... % kdir
            aap.tasklist.currenttask.settings.mask, ... % (mask)
            aap.tasklist.currenttask.settings.match, ... % (match)
            aap.tasklist.currenttask.settings.writeunwarpedEPI]; % (writeunwarpedEPI)
        
        % Fieldmap path
        FMdir = fullfile(aas_getsubjpath(aap, p), aap.directory_conventions.fieldmapsdirname);

        % This will work on all sessions (even those we have not selected)
        EPIdir = cell(size(aap.acq_details.sessions, 1));
        for s = aap.acq_details.selected_sessions
            % get files from stream
            EPIdir{s} = aas_getsesspath(aap,p,s);
        end
        
        % If we cannot find any images in the FMdir, check for subdirs and
        % move images
        FMfns = dir(fullfile(FMdir, '*.nii'));
        if isempty(FMfns)
            fldrDir = genpath(FMdir);
            while ~isempty(strtok(fldrDir, ':'))
                % Get each of the directories made by gendir
                [fldrCurr fldrDir] = strtok(fldrDir, ':');
                FMfns = dir(fullfile(fldrCurr, '*.nii'));
                for f = 1:length(FMfns)
                    unix(['mv ' fullfile(fldrCurr, FMfns(f).name) ...
                               ' ' FMdir]);
                end
            end
        end
        
        % @@@ CHANGE TO AUTOMATICALLY DETECT IMAGES USING aa4 tools...
        
        % @@@ THIS WILL NOT WORK FOR PROTOCOLS DIRECTLY CREATING THE VDM
        % @@@ ASK MARIJN KROES or ERNO HEMMANS
        % @@@ MEMORY GROUP!
        % @@@ PROBABLY NEEDS MASKING & SCALING! (-2048?) 0-2048, extrapolate, 
        % @@@ to phase value, & convert from siemens value to fieldmap val.
        % @@@ taking into account the TE!
        % @@@ PETER: look in header of sequence...
        %{
        <!--
            CHANGE ORDER THING TO CM, THEN CHANGE BACK
            TO ABSTRACT UNITS DEPNDING ON THE FoV
            start with 5 cm? ask Eelke
            FOV = resolution * by matrix size
            FOV ./ [order] = cm
            
            Also output jacobian map at the end, which should be similar to the fieldmap!
            
            
            cm -> order...
            
            ALTERNATIVE to get order:
            1) mask fieldmap data by BET
            2) take fourier transform of mask data = fieldmap in k space
            fftn
            collapse across 2 dimensions
            fftshift (or look at one side of data)
            3) look at the frequency peaks (maximal peak of spatial frequency of distortion)
            log scale...
            
            
            FIELDMAP:
            Look at the header information for possible scaling factor...
            
            COREGISTRATION:
            SCALE DATA BY 256 (min max)
            option...
            -->
        %}
        
        FieldMap_preprocess(FMdir,EPIdir,...
            pm_defs,...
            'session');
        
        outstream = {};
        
        % Rename VDM files to their correspondent run names
        if length(EPIdir) == 1
            VDMs = dir(fullfile(FMdir, 'vdm*.nii'));
            [~, fn, ext] = fileparts(VDMs.name);
            unix(['mv ' fullfile(FMdir, VDMs.name) ' ' ...
                fullfile(FMdir, [fn '_session1' ext])]);            
        end
        
        VDMs = dir(fullfile(FMdir, '*session*.nii'));

        for v = 1:length(VDMs)
            indx = strfind(VDMs(v).name, 'session');
            s = VDMs(v).name(indx+7:end-4); % Get number after 'session'
            s = str2double(s);
            
            % This gets the selected sessions!
            newfn = [VDMs(v).name(1:indx-1), ...
                aap.acq_details.sessions(aap.acq_details.selected_sessions(s)).name, '.nii'];
            unix(['mv ' fullfile(FMdir, VDMs(v).name)...
                ' ' fullfile(FMdir, newfn)]);
            outstream = [outstream fullfile(FMdir, newfn)];
        end

        if isempty(outstream)
            aas_log(aap, true, 'Could not find a fieldmap VDM after processing!')
        end
        
        aap=aas_desc_outputs(aap,p,'fieldmap',outstream);
        
        time_elapsed
end
