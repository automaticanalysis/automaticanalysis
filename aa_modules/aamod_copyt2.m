% AA module - converts T2 from DICOM and copies to central store
% If there is a T2 for this subject, copy it to central store
% Also makes 'structurals' directory
% Rhodri Cusack MRC CBU Cambridge Aug 2004

function [aap,resp]=aamod_copyt2(aap,task,i)

resp='';

switch task       
    case 'description'
        resp='T2 dicom to nifti and copying';
        
    case 'summary'
        if (length(aap.acq_details.subjects(i).t2)==0)
            resp=sprintf('No T2 for subject %s\n',aap.acq_details.subjects(i).mriname);
        else
            resp=sprintf('Converted T2 for subject %s \n', aas_getsubjname(aap,i));
        end;
        
    case 'report'
        
    case 'doit'
        
    [aap convertedfns dcmhdr]=aas_convertseries_fromstream(aap,i,'dicom_t2');
        
    % Save EXAMPLE dicom header (not all as previous code)
    subjpath=aas_getsubjpath(aap,i);

     % Save outputs?
    
    aap=aas_desc_outputs(aap,i,'t2',convertedfns);
            
    dcmhdrfn=fullfile(subjpath,'t2_dicom_header.mat');
    save(dcmhdrfn,'dcmhdr');
    aap=aas_desc_outputs(aap,i,'t2_dicom_header',dcmhdrfn);

 
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;