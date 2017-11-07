function [aap, resp]=aamod_diffusion_topup(aap,task,subjind,diffsessind)
resp='';

switch task
    case 'report'
    case 'doit'
        % Get inputs
        sliceaxes = {'ROW' 'COL'};
        pedirs = [+1 -1]; % flipping images to match Analyze
        nodif_topuptable=zeros(2,4);
        allfn='';
        for peind=1:2
            Yfn=aas_getfiles_bystream(aap,'diffusion_session_phaseencode_direction',[subjind diffsessind peind],'nodif');
            allfn=[allfn ' ' Yfn]; % used for fslmerge
            fn=aas_getfiles_bystream(aap,'diffusion_session_phaseencode_direction',[subjind diffsessind peind],'diffusion_dicom_header');
            DICOMHEADERS=load(fn); hdr = DICOMHEADERS.DICOMHEADERS{1};
            % PE
            nodif_topuptable(peind, cell_index(sliceaxes, deblank(hdr.InPlanePhaseEncodingDirection))) = ...
                pedirs(aas_get_numaris4_numval(hdr.CSAImageHeaderInfo,'PhaseEncodingDirectionPositive')+1);
            % Total Readout Time for FSL
            headerFields = fieldnames(hdr);
            nodif_topuptable(peind,4) = (hdr.(headerFields{strcmpi(headerFields,'NumberOfPhaseEncodingSteps')})-1)*hdr.echospacing;
        end
        
        % Output directory
        dsess=aas_getpath_bydomain(aap,'diffusion_session',[subjind diffsessind]);
        
        % Create merged file
        mergedfn=fullfile(dsess,'nodif_allpe.nii');
        cmd=['fslmerge -t ' mergedfn ' ' allfn];
        [s, w]=aas_runfslcommand(aap,cmd);       
        if s
            aas_log(aap,true,sprintf('Error running %s which was:\n%s',cmd,w));
        end
        
        % Write topup table
        apfn=fullfile(dsess,'acquisition_parameters.txt');
        fid=fopen(apfn,'w');
        for ind=1:size(nodif_topuptable,1)
            fprintf(fid,'%d %d %d %1.6f',nodif_topuptable(ind,:));
            fprintf(fid,'\n');
        end
        fclose(fid);
        
        % Now run topup
        outfn=fullfile(dsess,'topup_output');
        outfn_hifi_nodif=fullfile(dsess,'hifi_nodif');
        cmd=sprintf('topup --imain=%s --datain=%s --config=%s/etc/flirtsch/b02b0.cnf  --out=%s --iout=%s',mergedfn, apfn, aap.directory_conventions.fsldir, outfn, outfn_hifi_nodif);    
        aas_log(aap,false,sprintf('Running %s',cmd));        
        aas_runfslcommand(aap,cmd);
        
        % Describe outputs
        aap=aas_desc_outputs(aap,'diffusion_session',[subjind diffsessind],'topup_acquisition_parameters',apfn);
        aap=aas_desc_outputs(aap,'diffusion_session',[subjind diffsessind],'topup_output_movpar',fullfile(dsess,'topup_output_movpar.txt'));
        aap=aas_desc_outputs(aap,'diffusion_session',[subjind diffsessind],'topup_output_fieldcoef',spm_select('FPList',dsess, '^topup_output_fieldcoef.*'));
        % New corrected nodif
        aap=aas_desc_outputs(aap,'diffusion_session',[subjind diffsessind],'nodif',spm_select('FPList',dsess, '^hifi_nodif.*'));
end
end

