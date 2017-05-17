% AA module - Converts fieldmaps maps to NIFTI format
% Rhodri Cusack MRC CBU Cambridge Nov 2005

function [aap,resp]=aamod_convert_fieldmaps(aap,task,subj,sess)

resp='';

switch task
    case 'report'
    case 'doit'
        [aap, convertedfns, dcmhdr]=aas_convertseries_fromstream(aap,subj,sess,'dicom_fieldmap');
        
        sesspath=aas_getsesspath(aap,subj,sess);
        
        % Restructure outputs
        outstream = vertcat(convertedfns{:});
        aap=aas_desc_outputs(aap,subj,sess,'fieldmap',outstream);

        dcmhdr = dcmhdr{(numel(dcmhdr{2}) > numel(dcmhdr{1}))+1};
        dcmhdrfn=fullfile(sesspath,'fieldmap_dicom_header.mat');
        save(dcmhdrfn,'dcmhdr');
        aap=aas_desc_outputs(aap,subj,sess,'fieldmap_dicom_header',dcmhdrfn);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,true,sprintf('Unknown task %s',task));
end