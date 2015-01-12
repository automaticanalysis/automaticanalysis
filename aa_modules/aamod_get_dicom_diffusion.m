% This module finds all of the DICOM files associated with the diffusion
% using aas_listdicomfiles, and copies them into the session directory of
% this module, either across the local filesystem or from s3. It then
% creates the output stream.
% function aap=aamod_get_dicom_structural(aap,task,subjind,diffsessind)

function [aap resp]=aamod_get_dicom_diffusion(aap,task,subjind,diffsessind)
global aaworker
resp='';

switch task
    case 'report'
    case 'doit'
        dsesspth= aas_getpath_bydomain(aap,'diffusion_session',[subjind diffsessind]);
        [aap fns]=aas_listdicomfiles(aap,subjind,aap.acq_details.subjects(subjind).diffusion_seriesnumbers(diffsessind));
        
        outfns={};
        for fnind=1:length(fns)
            copyfile(fns{fnind},dsesspth);
            [pth nme ext]=fileparts(fns{fnind});
            outfns{end+1}=[nme ext];
        end;
        
        aap=aas_desc_outputs(aap,'diffusion_session',[subjind diffsessind],'dicom_diffusion',outfns);
end
end

