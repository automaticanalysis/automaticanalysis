% AA module
% Convert the fieldmap images (2 mag and 1 phase) into a Voxel Displacement
% Map (VDM) using the FieldMap toolbox of SPM

function [aap,resp]=aamod_fieldmap2VDM(aap,task,subj)

resp='';

switch task
    case 'domain'
        resp='subject';  % this module needs to be run once per subject
        
    case 'description'
        resp='SPM5 align';
        
    case 'summary'
        subjpath=aas_getsubjpath(subj);
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

        % retrieve actual values (based on the first EPI) --> TODO: session-specific fieldmaps
        EPI_DICOMHEADERS = load(aas_getimages_bystream(aap,subj,1,'epi_dicom_header')); EPI_DICOMHEADERS = EPI_DICOMHEADERS.DICOMHEADERS{1};
        if isfield(EPI_DICOMHEADERS,'NumberOfPhaseEncodingSteps') && isfield(EPI_DICOMHEADERS,'echospacing')
            pm_defs(4) = EPI_DICOMHEADERS.NumberOfPhaseEncodingSteps*EPI_DICOMHEADERS.echospacing*1000;
        else
            aas_log(aap,true,'ERROR:Field for number of phase encoding steps and/or echospacing not found!\nERROR:You may need to rerun aamod_convert_epis.');
        end
        try % te1 and te2
            FM_DICOMHEADERS=load(aas_getfiles_bystream(aap,subj,'fieldmap_dicom_header'));
            dcmhdrs = cell2mat(FM_DICOMHEADERS.dcmhdr);
            tes = sort(unique([dcmhdrs.EchoTime]),'ascend');
            if numel(tes) > 2, error('wrong!'); end
            pm_defs(1:2) = tes;
        catch
            aas_log(aap,0,'Error during retrieving Fieldmap Echo Times - Defaults are used!');
        end
        if ~all(pm_defs([1 2 4])) % retrieval error and no defaults specified
            aas_log(aap,1,'Error during retrieving Timings and no Default Timings found!');
        end
        
        % Fieldmap path
        FMdir = fullfile(aas_getsubjpath(aap, subj), aap.directory_conventions.fieldmapsdirname);
        
        % The folder fieldmaps must exist...
        if ~exist(FMdir, 'dir')
            mkdir(FMdir)
        else
            % remove previous vdms
            delete(fullfile(FMdir, 'vdm*.nii'));
        end
        
        FMfn = aas_getfiles_bystream(aap,subj,'fieldmap');
        
        % This will work on all sessions (even those we have not selected)
        EPIdir = cell(size(aap.acq_details.sessions, 1));
        for s = aap.acq_details.selected_sessions
            % get files from stream
            EPIdir{s} = aas_getsesspath(aap,subj,s);
        end
        
        % If we cannot find any images in the FMdir, move images there...
        FMfns = dir(fullfile(FMdir, '*.nii'));
        if isempty(FMfns)
            for f = 1:size(FMfn,1)
                aas_shell(['mv ' squeeze(FMfn(f,:)) ' ' FMdir]);
            end
        end
        FieldMap_preprocess(FMdir,EPIdir,...
            pm_defs,...
            'session');
        
        outstream = {};
        
        % Rename VDM files to their correspondent run names
        if length(EPIdir) == 1
            VDMs = dir(fullfile(FMdir, 'vdm*.nii'));
            [junk, fn, ext] = fileparts(VDMs.name);
            aas_shell(['mv ' fullfile(FMdir, VDMs.name) ' ' ...
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
            aas_shell(['mv ' fullfile(FMdir, VDMs(v).name)...
                ' ' fullfile(FMdir, newfn)]);
            outstream = [outstream fullfile(FMdir, newfn)];
        end
        
        if isempty(outstream)
            aas_log(aap, true, 'Could not find a fieldmap VDM after processing!')
        end
        
        aap=aas_desc_outputs(aap,subj,'fieldmap',outstream);
        
end
