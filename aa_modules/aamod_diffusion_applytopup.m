function [aap, resp]=aamod_diffusion_applytopup(aap,task,subjind,diffsessind)
resp='';

switch task
    case 'report'
    case 'doit'
        
        % Get nii filenames from stream
        allfns='';
        for peind=1:2
            diffinput=aas_getfiles_bystream(aap,'diffusion_session_phaseencode_direction',[subjind diffsessind peind],'diffusion_data');
            allfns=[allfns ',' diffinput];
        end;
        allfns(1) = '';
        % Output directory
        dsess=aas_getpath_bydomain(aap,'diffusion_session',[subjind diffsessind]);
         
        acqparms=aas_getfiles_bystream(aap,'diffusion_session',[subjind diffsessind],'topup_acquisition_parameters');
        topup_output_base=char(aas_getfiles_bystream(aap,'diffusion_session',[subjind diffsessind],'topup_output_fieldcoef'),...
            aas_getfiles_bystream(aap,'diffusion_session',[subjind diffsessind],'topup_output_movpar'));
        topup_output_base = topup_output_base(1,topup_output_base(1,:) == topup_output_base(2,:)); topup_output_base(end) = '';
                
        % Now apply topup
        cmd = sprintf('applytopup --imain=%s  --datain=%s  --inindex=1,2 --topup=%s --out=%s',allfns,acqparms,topup_output_base,fullfile(dsess,'topup_dwi'));
        aas_log(aap,false,sprintf('Running %s',cmd));        
        aas_runfslcommand(aap,cmd);
        
        % Describe outputs
        outfn = spm_select('FPList',dsess,'^topup_dwi.*');
        % clean data
        V = spm_vol(outfn);
        Y = spm_read_vols(V);
        Y(Y<0) = 0;
        for vol = 1:numel(V)
            V(vol).fname = spm_file(V(vol).fname,'suffix',sprintf('_%03d',V(vol).n(1)));
            Vout(vol) = spm_write_vol(V(vol),Y(:,:,:,vol));
        end
        delete(outfn);
        if aap.options.NIFTI4D
            spm_file_merge(Vout,outfn);
            fname = {Vout.fname};
            delete(fname{:});
            Vout = Vout(1); Vout.fname = outfn;
        end
        % data
        aap=aas_desc_outputs(aap,'diffusion_session',[subjind diffsessind],'diffusion_data',{Vout.fname});
        % save bvals bvecs from 1st
        aap=aas_desc_outputs(aap,'diffusion_session',[subjind diffsessind],'bvals',...
            aas_getfiles_bystream(aap,'diffusion_session_phaseencode_direction',[subjind diffsessind 1],'bvals'));
        aap=aas_desc_outputs(aap,'diffusion_session',[subjind diffsessind],'bvecs',...
            aas_getfiles_bystream(aap,'diffusion_session_phaseencode_direction',[subjind diffsessind 1],'bvecs'));
end
end

