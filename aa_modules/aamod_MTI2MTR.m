function [aap,resp]=aamod_MTI2MTR(aap,task,subj,sess)

resp='';

switch task
    case 'report'
        
    case 'doit'
        %% Init
        global defaults
        flags = defaults.coreg;
        flags.write.interp = 1;
        flags.write.which = [1 0]; % do not (re)write target and mean
        
        %% Data
        % MTI path
        mtifn = aas_getfiles_bystream(aap,'special_session',[subj,sess],'MTI');
        mtihdr = load(aas_getfiles_bystream(aap,'special_session',[subj,sess],'MTI_dicom_header'));
        
        % 2 series expected: baseline and MT
        mtihdr = {mtihdr.dcmhdr{1}.SeriesDescription mtihdr.dcmhdr{2}.SeriesDescription};
        baseindex = cell_index(mtihdr,'baseline');
        mtibaseline = deblank(mtifn(baseindex,:));
        mtiMT = deblank(mtifn(~(baseindex-1)+1,:));
        
        %% Coregister
        % Coregister MT to baseline
        x = spm_coreg(spm_vol(mtibaseline), spm_vol(mtiMT), flags.estimate);
        % Set the new space for MT
        spm_get_space(mtiMT, spm_matrix(x)\spm_get_space(mtiMT));
                
        aas_log(aap,false,sprintf(['\tMTI to baseline realignment parameters:\n' ...
            '\tx: %0.4f   y: %0.4f   z: %0.4f   p: %0.4f   r: %0.4f   j: %0.4f'], ...
            x(1), x(2), x(3), x(4), x(5), x(6)))

        % Reslice MT
        spm_reslice(strvcat(mtibaseline,mtiMT),flags.write);
        [pth, fname, ext] = fileparts(mtiMT);
        mtiMT = fullfile(pth, [flags.write.prefix fname ext]);
        
        %% Skull-strip MT and baseline
        bmtibaseline = spm_file(mtibaseline,'prefix','b');
        aas_runfslcommand(aap,sprintf('bet %s %s -f 0.6',mtibaseline,bmtibaseline));
        bmtiMT = spm_file(mtiMT,'prefix','b');
        aas_runfslcommand(aap,sprintf('bet %s %s -f 0.6',mtiMT,bmtiMT));
        
        %% Compute MTR
        R1 = spm_read_vols(spm_vol(bmtiMT));
        R2 = spm_read_vols(spm_vol(bmtibaseline));
        MTR = ((R2 - R1)./R2) * 100;    % MTR calculation
        MTR(MTR < 0) = 0;               % thresholds the image to remove negative voxels
        MTR(MTR > 90) = 0;              % thresholds the image to remove positive outliers
        MTR(isnan(MTR)) = 0;

        V = spm_vol(bmtiMT);
        mtrfn = fullfile(aas_getsesspath(aap,subj,sess),'mti_MTR.nii');
        nifti_write(mtrfn,MTR,'MTR',V);        
        aap=aas_desc_outputs(aap,'special_session',[subj,sess],'MTR',mtrfn);
        aap=aas_desc_outputs(aap,'special_session',[subj,sess],'MTI_baseline',mtibaseline);
end