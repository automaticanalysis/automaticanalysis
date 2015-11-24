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
        
        % Defaults
        pm_defaults;
        pm_def.sessname='session';
        pm_def.EPI_BASED_FIELDMAPS = aas_getsetting(aap,'epifm');        
        pm_def.K_SPACE_TRAVERSAL_BLIP_DIR = aas_getsetting(aap,'kdir');
        if pm_def.EPI_BASED_FIELDMAPS==1
            pm_def.INPUT_DATA_FORMAT='RI';
            pm_def.MASKBRAIN=0;
        elseif pm_def.EPI_BASED_FIELDMAPS==0
            pm_def.INPUT_DATA_FORMAT='PM';
            pm_def.MASKBRAIN=aas_getsetting(aap,'mask');
        end
        pm_def.match_vdm = aas_getsetting(aap,'match');
        pm_def.write_unwarped = aas_getsetting(aap,'writeunwarpedEPI');
        
        % retrieve actual values (based on the first EPI) --> TODO: session-specific fieldmaps
        EPI_DICOMHEADERS = load(aas_getimages_bystream(aap,subj,1,'epi_dicom_header')); EPI_DICOMHEADERS = EPI_DICOMHEADERS.DICOMHEADERS{1};
        if isfield(EPI_DICOMHEADERS,'NumberOfPhaseEncodingSteps') && isfield(EPI_DICOMHEADERS,'echospacing')
            pm_def.TOTAL_EPI_READOUT_TIME = EPI_DICOMHEADERS.NumberOfPhaseEncodingSteps*EPI_DICOMHEADERS.echospacing*1000;
        else
            aas_log(aap,true,'ERROR:Field for number of phase encoding steps and/or echospacing not found!\nERROR:You may need to rerun aamod_convert_epis.');
        end
        
        try % te1 and te2
            FM_DICOMHEADERS=load(aas_getfiles_bystream(aap,subj,'fieldmap_dicom_header'));
            dcmhdrs = cell2mat(FM_DICOMHEADERS.dcmhdr);
            tes = sort(unique([dcmhdrs.EchoTime]),'ascend');
            if numel(tes) ~= 2, error('Inappropriate fieldmap header!'); end
            pm_def.SHORT_ECHO_TIME = tes(1);
            pm_def.LONG_ECHO_TIME = tes(2);
        catch E
            pm_def.SHORT_ECHO_TIME = aas_getsetting(aap,'te1');
            pm_def.LONG_ECHO_TIME = aas_getsetting(aap,'te2');
            aas_log(aap,false,sprintf('WARNING: Error during retrieving Fieldmap Echo Times: %s\nWARNING: Defaults are used!',E.message));
        end
        
        % Fieldmap path
        FMdir = fullfile(aas_getsubjpath(aap, subj), aap.directory_conventions.fieldmapsdirname);
%         delete(fullfile(FMdir, '*.nii')); % remove previous vdms
        FM = aas_getfiles_bystream(aap,subj,'fieldmap');
        for f = 1:size(FM,1)
            aas_shell(['cp ' squeeze(FM(f,:)) ' ' FMdir]);
            FMfn{f} = spm_file(FM(f,:),'path',FMdir);
        end
        scphase=FieldMap('Scale',FMfn{3});
        fm_imgs=char(scphase.fname,FMfn{1});
        
        % This will work on all sessions (even those we have not selected)
        epi_img = cell(size(aap.acq_details.sessions, 1));
        for s = aap.acq_details.selected_sessions
            % get files from stream
            epi_img{s} = spm_file(aas_getfiles_bystream(aap,subj,s,'epi'),'number',1);
        end
        
%         FieldMap_create(fm_imgs,epi_img,pm_def);
        
        % Rename VDM files to their correspondent run names
        if numel(epi_img) == 1
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
