% AA module
% Convert the fieldmap images (2 mag and 1 phase) into a Voxel Displacement
% Map (VDM) using the FieldMap toolbox of SPM

function [aap,resp]=aamod_fieldmap2VDM(aap,task,subj,sess)

resp='';

switch task
    case 'report'
%         if aas_getsetting(aap,'writeunwarpedEPI',sess)
%             spm_file(aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,[subj,sess],aas_getstreams(aap,'input',1)),'prefix','u')
%         end
    case 'doit'
        domain = aap.tasklist.currenttask.domain;
        
        % Fieldmap path
        FMdir = fullfile(aas_getsesspath(aap, subj,sess), aap.directory_conventions.fieldmapsdirname);
        delete(fullfile(FMdir, '*.nii')); % remove previous vdms
        
        % Fieldmaps
        FM = aas_getfiles_bystream(aap,domain,[subj,sess],'fieldmap');
        for f = 1:size(FM,1)
            aas_shell(['cp ' squeeze(FM(f,:)) ' ' FMdir]);
            FMfn{f} = spm_file(FM(f,:),'path',FMdir);
        end
        switch size(FM,1)
            case 8 % 2* (mag & phase & real & imag) (GE)
                job.data.realimag.shortreal = FMfn(3);
                job.data.realimag.shortimag = FMfn(4);
                job.data.realimag.longreal = FMfn(7);
                job.data.realimag.longimag = FMfn(8);
            case 4 % (real & real & imag & imag) (Philips)
                job.data.realimag.shortreal = FMfn(1);
                job.data.realimag.shortimag = FMfn(3);
                job.data.realimag.longreal = FMfn(2);
                job.data.realimag.longimag = FMfn(4);
            case 3 % presubphase + 2*magnitude (Siemens)
                job.data.presubphasemag.phase = FMfn(3);
                job.data.presubphasemag.magnitude = FMfn(1:2);
            case 1 % precalcfieldmap
                job.data.precalcfieldmap.precalcfieldmap = FMfn;
        end
        
        % EPI
        job.session.epi{1} = aas_getfiles_bystream(aap,domain,[subj,sess],aas_getstreams(aap,'input',1));
        
        % Flags
        job.defaults.defaultsval = FieldMap('SetParams');
        job.defaults.defaultsval.et = cell2mat(job.defaults.defaultsval.et);
        job.defaults.defaultsval.mflags.template = cellstr(job.defaults.defaultsval.mflags.template);
        job.sessname = 'session';
        job.anat = [];
        
        job.defaults.defaultsval.epifm = aas_getsetting(aap,'epifm',sess);
        job.defaults.defaultsval.blipdir = aas_getsetting(aap,'kdir',sess);
        job.defaults.defaultsval.maskbrain = aas_getsetting(aap,'mask',sess);
        job.matchvdm = aas_getsetting(aap,'match',sess);
        job.writeunwarped = aas_getsetting(aap,'writeunwarpedEPI',sess);
        
        % EPI TotalEPIReadoutTime (based on the first EPI)
        EPI_DICOMHEADERS = load(aas_getfiles_bystream(aap,domain,[subj,sess],aas_getstreams(aap,'input',2))); EPI_DICOMHEADERS = EPI_DICOMHEADERS.DICOMHEADERS{1};
        if isfield(EPI_DICOMHEADERS,'NumberofPhaseEncodingSteps') && isfield(EPI_DICOMHEADERS,'echospacing') % <=SPM8
            job.defaults.defaultsval.tert = EPI_DICOMHEADERS.NumberofPhaseEncodingSteps*EPI_DICOMHEADERS.echospacing*1000;
        elseif isfield(EPI_DICOMHEADERS,'NumberOfPhaseEncodingSteps') && isfield(EPI_DICOMHEADERS,'echospacing') % >=SPM12
            job.defaults.defaultsval.tert = EPI_DICOMHEADERS.NumberOfPhaseEncodingSteps*EPI_DICOMHEADERS.echospacing*1000;
        else
            job.defaults.defaultsval.tert = aas_getsetting(aap,'tert',sess);
            aas_log(aap,false,sprintf('WARNING: Field for number of phase encoding steps and/or echospacing not found!\nWARNING: Manual setting is used: %1.3f ms.',job.defaults.defaultsval.tert));
            if isempty(job.defaults.defaultsval.tert), aas_log(aap,true,'ERROR: No value is specified'); end
        end
        
        if ~isfield(job.data,'precalcfieldmap')
            try % Fieldmap EchoTimes
                FM_DICOMHEADERS=load(aas_getfiles_bystream(aap,domain,[subj,sess],'fieldmap_dicom_header'));
                if iscell(FM_DICOMHEADERS.dcmhdr)
                    assert(numel(FM_DICOMHEADERS.dcmhdr{1,1})==1, 'unexpected input');
                    tes = [FM_DICOMHEADERS.dcmhdr{1,1}.EchoTime1*1000,FM_DICOMHEADERS.dcmhdr{1,1}.EchoTime2*1000]; % in ms
                else
                    tes = sort(unique(cellfun(@(x) x.volumeTE,FM_DICOMHEADERS.dcmhdr)),'ascend')*1000; % in ms
                end
                if numel(tes) ~= 2, tes = sort(unique(cellfun(@(x) x.EchoTime,FM_DICOMHEADERS.dcmhdr)),'ascend'); end % try original (backward compatibility)
                if numel(tes) ~= 2, aas_log(aap,true,'Inappropriate fieldmap header!'); end
                job.defaults.defaultsval.et = tes;
            catch E
                job.defaults.defaultsval.et = [aas_getsetting(aap,'te1',sess) aas_getsetting(aap,'te2',sess)];
                aas_log(aap,false,sprintf('WARNING: Error during retrieving Fieldmap Echo Times: %s\nWARNING: Manual settings are used: %1.3f ms and %1.3f ms!',E.message,job.defaults.defaultsval.et));
                if isempty(job.defaults.defaultsval.et), aas_log(aap,true,'ERROR: No value is specified'); end
            end
        end
        
        FieldMap_Run(job);
        
        % Rename VDM files to their correspondent run names
        
        VDM = spm_select('FPList',FMdir,'^vdm.*nii');
        if isempty(VDM)
            aas_log(aap, true, 'ERROR: Could not find a fieldmap VDM after processing!')
        end
        
        outstream = spm_file(VDM,'suffix',['_' aap.acq_details.sessions(sess).name]);
        movefile(VDM,outstream);
        aap=aas_desc_outputs(aap,domain,[subj,sess],'fieldmap',outstream);                
end
