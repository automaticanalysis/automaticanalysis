% Extract the reference(s) image(s) (T2 image with b-value of 0), called
% nodif

function [aap resp]=aamod_diffusion_topup(aap,task,subjind,diffsessind)
resp='';

switch task
    case 'report'
    case 'doit'
        % Get inputs
        nodif_allpe=[];
        allfn='';
        for peind=1:2
            Yfn=aas_getfiles_bystream(aap,'diffusion_session_phaseencode_direction',[subjind diffsessind peind],'nodif');
            allfn=[allfn ' ' Yfn]; % used for fslmerge
% In future, could include checks of topup table using DICOM headers
%            fn=aas_getfiles_bystream(aap,'diffusion_session_phaseencode_direction',[subjind diffsessind peind],'diffusion_dicom_header');
%            nodif_header(peind)=load(fn);
        end;
        
        % Output directory
        dsess=aas_getpath_bydomain(aap,'diffusion_session',[subjind diffsessind]);
        
        % Create merged file
        mergedfn=fullfile(dsess,'nodif_allpe.nii');
        cmd=['fslmerge -t ' mergedfn ' ' allfn];
        [s w]=aas_runfslcommand(aap,cmd);       
        if s
            aas_log(aap,true,sprintf('Error running %s which was:\n%s',cmd,w));
        end;
        
        % Write topup table
        apfn=fullfile(dsess,'acquisition_parameters.txt');
        fid=fopen(apfn,'w');
        for ind=1:length(aap.tasklist.currenttask.settings.topuptable.topuprow)
            fprintf(fid,'%f ',aap.tasklist.currenttask.settings.topuptable.topuprow{ind});
            fprintf(fid,'\n');
        end;
        fclose(fid);
        
        % Now run topup
        outfn=fullfile(dsess,'topup_output');
        outfn_hifi_nodif=fullfile(dsess,'hifi_nodif');
        cmd=sprintf('topup --imain=%s --datain=%s --config=b02b0.cnf  --out=%s --iout=%s',mergedfn, apfn,outfn,outfn_hifi_nodif);    
        aas_log(aap,false,sprintf('Running %s',cmd));        
        [s w]=aas_runfslcommand(aap,cmd);
        if s
            aas_log(aap,true,sprintf('Error %s',w));
        end;
        
        % Check format of output
        fieldcoef_fn=dir([outfn '_fieldcoef.ni*']);
        hifi_nodif_fn=dir([outfn_hifi_nodif '.ni*']);
        
        % Describe outputs
        aap=aas_desc_outputs(aap,'diffusion_session',[subjind diffsessind],'topup_acquisition_parameters',apfn);
        aap=aas_desc_outputs(aap,'diffusion_session',[subjind diffsessind],'topup_output_movpar',fullfile(dsess,'topup_output_movpar.txt'));
        aap=aas_desc_outputs(aap,'diffusion_session',[subjind diffsessind],'topup_output_fieldcoef',fullfile(dsess,fieldcoef_fn(1).name));
        % New corrected nodif
        aap=aas_desc_outputs(aap,'diffusion_session',[subjind diffsessind],'nodif',fullfile(dsess,hifi_nodif_fn(1).name));
end
end

