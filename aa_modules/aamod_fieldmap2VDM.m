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
        
        % The folder fieldmaps must exist...
        if ~exist(FMdir, 'dir')
            mkdir(FMdir)
        end
        
        FMfn = aas_getfiles_bystream(aap,p,'fieldmap');
        
        % This will work on all sessions (even those we have not selected)
        EPIdir = cell(size(aap.acq_details.sessions, 1));
        for s = aap.acq_details.selected_sessions
            % get files from stream
            EPIdir{s} = aas_getsesspath(aap,p,s);
        end
        
        % If we cannot find any images in the FMdir, move images there...
        FMfns = dir(fullfile(FMdir, '*.nii'));
        if isempty(FMfns)
            for f = 1:size(FMfn,1)
                unix(['mv ' squeeze(FMfn(f,:)) ' ' FMdir])
            end
        end
        
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
                aap.acq_details.sessions(aap.acq_details.selected_sessions(s)).name,'.nii'];
            unix(['mv ' fullfile(FMdir, VDMs(v).name)...
                ' ' fullfile(FMdir, newfn)]);
            outstream = [outstream fullfile(FMdir, newfn)];
        end
        
        if isempty(outstream)
            aas_log(aap, true, 'Could not find a fieldmap VDM after processing!')
        end
        
        aap=aas_desc_outputs(aap,p,'fieldmap',outstream);
        
end