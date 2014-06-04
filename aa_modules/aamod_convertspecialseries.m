% AA module - Converts special series to NIFTI format
% Rhodri Cusack MRC CBU Cambridge Nov 2005

function [aap,resp]=aamod_convertspecialseries(aap,task,subjI)

resp='';

switch task
    case 'domain'
        resp='subject';   % this module needs to be run once per subject

    case 'description'
        resp='Convert special series';

    case 'summary'
        somethingtoconvert=false;
        if (length(aap.acq_details.specialseries)~=0)
            if (length(aap.acq_details.specialseries{subjI})==0)
                somethingtoconvert=true;
            end;
        end;
        if (somethingtoconvert)
            resp=sprintf('No special series for subject %s\n',aap.acq_details.subjects(subjI).mriname);
        else
            resp=sprintf('Converted special series for subject %s to %s\n', aap.acq_details.subjects(subjI).mriname,aap.directory_conventions.centralstore_structurals);
        end;

    case 'report'
    case 'doit'
        
        [aap convertedfns dcmhdr] = aas_convertseries_fromstream(aap, subjI, 'dicom_specialseries'); 
        
        
        aap = aas_desc_outputs(aap, subjI, 'specialseries', convertedfns);
        
    case 'checkrequirements'

    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
