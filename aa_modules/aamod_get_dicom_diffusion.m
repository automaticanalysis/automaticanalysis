% This module finds all of the DICOM files associated with the diffusion
% using aas_listdicomfiles, and copies them into the session directory of
% this module, either across the local filesystem or from s3. It then
% creates the output stream.
% function aap=aamod_get_dicom_structural(aap,task,subjind,diffsessind)

function [aap resp]=aamod_get_dicom_diffusion(aap,task,subj,sess)
global aaworker
resp='';

switch task
    case 'report'
    case 'doit'
        dsesspth= aas_getpath_bydomain(aap,'diffusion_session',[subj sess]);
        [d, mriser] = aas_get_series(aap,'diffusion',subj,sess);
        
        [aap fns]=aas_listdicomfiles(aap,[subj d],mriser);
        
        outfns={};
        for fnind=1:length(fns)
            copyfile(fns{fnind},dsesspth);
            [pth nme ext]=fileparts(fns{fnind});
            outfns{end+1}=[nme ext];
        end;
        
        aap=aas_desc_outputs(aap,'diffusion_session',[subj sess],'dicom_diffusion',outfns);
end
end

