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
        % Fieldmap path
        FMdir = fullfile(aas_getsubjpath(aap, subj), aap.directory_conventions.fieldmapsdirname);
        delete(fullfile(FMdir, '*.nii')); % remove previous vdms
        
        % Flags
        job.defaults.defaultsval = FieldMap('SetParams');
        job.defaults.defaultsval.mflags.template = cellstr(job.defaults.defaultsval.mflags.template);
        job.sessname = 'session';
        job.anat = [];
        
        job.defaults.defaultsval.epifm = aas_getsetting(aap,'epifm');
        job.defaults.defaultsval.blipdir = aas_getsetting(aap,'kdir');
        job.defaults.defaultsval.maskbrain = aas_getsetting(aap,'mask');
        job.matchvdm = aas_getsetting(aap,'match');
        job.writeunwarped = aas_getsetting(aap,'writeunwarpedEPI');
        
        % EPI TotalEPIReadoutTime (based on the first EPI) --> TODO: session-specific fieldmaps
        EPI_DICOMHEADERS = load(aas_getimages_bystream(aap,subj,1,'epi_dicom_header')); EPI_DICOMHEADERS = EPI_DICOMHEADERS.DICOMHEADERS{1};
        if isfield(EPI_DICOMHEADERS,'NumberOfPhaseEncodingSteps') && isfield(EPI_DICOMHEADERS,'echospacing')
            job.defaults.defaultsval.tert = EPI_DICOMHEADERS.NumberOfPhaseEncodingSteps*EPI_DICOMHEADERS.echospacing*1000;
        else
            aas_log(aap,true,'ERROR:Field for number of phase encoding steps and/or echospacing not found!\nERROR:You may need to rerun aamod_convert_epis.');
        end
        
        try % Fieldmap EchoTimes
            FM_DICOMHEADERS=load(aas_getfiles_bystream(aap,subj,'fieldmap_dicom_header'));
            dcmhdrs = cell2mat(FM_DICOMHEADERS.dcmhdr);
            tes = sort(unique([dcmhdrs.EchoTime]),'ascend');
            if numel(tes) ~= 2, error('Inappropriate fieldmap header!'); end
            job.defaults.defaultsval.et = tes;
        catch E
            job.defaults.defaultsval.et = [aas_getsetting(aap,'te1') aas_getsetting(aap,'te2')];
            aas_log(aap,false,sprintf('WARNING: Error during retrieving Fieldmap Echo Times: %s\nWARNING: Defaults are used!',E.message));
        end
        
        % Fieldmaps
        FM = aas_getfiles_bystream(aap,subj,'fieldmap');
        for f = 1:size(FM,1)
            aas_shell(['cp ' squeeze(FM(f,:)) ' ' FMdir]);
            FMfn{f} = spm_file(FM(f,:),'path',FMdir);
        end
        job.data.presubphasemag.phase = FMfn(3);
        job.data.presubphasemag.magnitude = FMfn(1:2);
        
        % EPI: This will work on all sessions (even those we have not selected)
        for s = aap.acq_details.selected_sessions
            % get files from stream
            job.session(s).epi{1} = spm_file(aas_getfiles_bystream(aap,subj,s,'epi'),'number',1);
        end
        
        FieldMap_Run(job);
        
        % Rename VDM files to their correspondent run names
        if numel(job.session) == 1
            VDMs = dir(fullfile(FMdir, 'vdm*.nii'));
            [junk, fn, ext] = fileparts(VDMs.name);
            aas_shell(['mv ' fullfile(FMdir, VDMs.name) ' ' ...
                fullfile(FMdir, [fn '_session1' ext])]);
        end
        
        VDMs = dir(fullfile(FMdir, '*session*.nii'));
        
        outstream = {};
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
