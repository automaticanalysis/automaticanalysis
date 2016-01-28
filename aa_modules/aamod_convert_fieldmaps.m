% AA module - Converts fieldmaps maps to NIFTI format
% Rhodri Cusack MRC CBU Cambridge Nov 2005

function [aap,resp]=aamod_convert_fieldmaps(aap,task,subj)

resp='';

switch task       
    case 'description'
        resp='Fieldmaps dicom to nifti and copying';
        

    case 'report'
    case 'doit'
        
    [aap convertedfns dcmhdr]=aas_convertseries_fromstream(aap,subj,'dicom_fieldmap');
    
    dcmhdr = dcmhdr{(numel(dcmhdr{2}) > numel(dcmhdr{1}))+1};
    
    subjpath=aas_getsubjpath(aap,subj);

    outstream = {};
    % Restructure outputs!
    for c = 1:length(convertedfns)
        outstream = [outstream; convertedfns{c}];
    end
    
    % Save outputs?
    aap=aas_desc_outputs(aap,subj,'fieldmap',outstream);
            
    dcmhdrfn=fullfile(subjpath,'fieldmap_dicom_header.mat');
    save(dcmhdrfn,'dcmhdr');
    aap=aas_desc_outputs(aap,subj,'fieldmap_dicom_header',dcmhdrfn);
    
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end