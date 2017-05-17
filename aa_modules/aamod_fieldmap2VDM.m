% AA module
% Convert the fieldmap images (2 mag and 1 phase) into a Voxel Displacement
% Map (VDM) using the FieldMap toolbox of SPM

function [aap,resp]=aamod_fieldmap2VDM(aap,task,subj,sess)

resp='';

switch task
    case 'report'
        
    case 'doit'
        % Fieldmap path
        FMdir = fullfile(aas_getsesspath(aap, subj,sess), aap.directory_conventions.fieldmapsdirname);
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
        
        % EPI TotalEPIReadoutTime (based on the first EPI)
        EPI_DICOMHEADERS = load(aas_getimages_bystream(aap,subj,sess,'epi_dicom_header')); EPI_DICOMHEADERS = EPI_DICOMHEADERS.DICOMHEADERS{1};
        if isfield(EPI_DICOMHEADERS,'NumberofPhaseEncodingSteps') && isfield(EPI_DICOMHEADERS,'echospacing') % <=SPM8
            job.defaults.defaultsval.tert = EPI_DICOMHEADERS.NumberofPhaseEncodingSteps*EPI_DICOMHEADERS.echospacing*1000;
        elseif isfield(EPI_DICOMHEADERS,'NumberOfPhaseEncodingSteps') && isfield(EPI_DICOMHEADERS,'echospacing') % >=SPM12
            job.defaults.defaultsval.tert = EPI_DICOMHEADERS.NumberOfPhaseEncodingSteps*EPI_DICOMHEADERS.echospacing*1000;
        else
            aas_log(aap,true,'ERROR:Field for number of phase encoding steps and/or echospacing not found!\nERROR:You may need to rerun aamod_convert_epis.');
        end
        
        try % Fieldmap EchoTimes
            FM_DICOMHEADERS=load(aas_getfiles_bystream(aap,subj,sess,'fieldmap_dicom_header'));
            dcmhdrs = cell2mat(FM_DICOMHEADERS.dcmhdr);
            tes = sort(unique([dcmhdrs.volumeTE]),'ascend')*1000; % in ms
            if numel(tes) ~= 2, tes = sort(unique([dcmhdrs.EchoTime]),'ascend'); end % try original (backward compatibility)
            if numel(tes) ~= 2, aas_log(aap,true,'Inappropriate fieldmap header!'); end
            job.defaults.defaultsval.et = tes;
        catch E
            job.defaults.defaultsval.et = [aas_getsetting(aap,'te1') aas_getsetting(aap,'te2')];
            aas_log(aap,false,sprintf('WARNING: Error during retrieving Fieldmap Echo Times: %s\nWARNING: Defaults are used!',E.message));
        end
        
        % Fieldmaps
        FM = aas_getfiles_bystream(aap,subj,sess,'fieldmap');
        for f = 1:size(FM,1)
            aas_shell(['cp ' squeeze(FM(f,:)) ' ' FMdir]);
            FMfn{f} = spm_file(FM(f,:),'path',FMdir);
        end
        switch size(FM,1)
            case 3 % presubphase + 2*magnitude
                job.data.presubphasemag.phase = FMfn(3);
                job.data.presubphasemag.magnitude = FMfn(1:2);
            case 1 % precalcfieldmap
                job.data.precalcfieldmap.precalcfieldmap = FM(fn);
        end
        
        % EPI
        job.session.epi{1} = spm_file(aas_getfiles_bystream(aap,subj,sess,'epi'),'number',1);
        
        FieldMap_Run(job);
        
        % Rename VDM files to their correspondent run names
        
        VDM = spm_select('FPList',FMdir,'^vdm.*nii');
        if isempty(VDM)
            aas_log(aap, true, 'ERROR: Could not find a fieldmap VDM after processing!')
        end
        
        outstream = spm_file(VDM,'suffix',['_' aap.acq_details.sessions(sess).name]);
        movefile(VDM,outstream);
        aap=aas_desc_outputs(aap,subj,sess,'fieldmap',outstream);        
end
