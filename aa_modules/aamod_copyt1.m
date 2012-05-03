% AA module - converts T1 structural from DICOM and copies to central store
% Also makes 'structurals' directory
% Original version: Rhodri Cusack MRC CBU Cambridge Aug 2004

function [aap,resp]=aamod_copyt1(aap,task,i)

resp='';

switch task
    case 'description'
        resp='T1 dicom to nifti and copying';

    case 'summary'
        if (length(aap.acq_details.subjects(i).structural)==0)
            resp = sprintf('No T1 structural image for subject %s\n',aap.acq_details.subjects(i).mriname);
        else
            resp = sprintf('Converted T1 structural image for subject %s \n', aas_getsubjname(aap,i));
        end;

    case 'report'
    case 'doit'

    [aap convertedfns dcmhdr] = aas_convertseries_fromstream(aap,i,'dicom_t1');

    % Save EXAMPLE dicom header (not all as previous code)
    subjpath = aas_getsubjpath(aap,i);

     % Save outputs?

    aap=aas_desc_outputs(aap,i,'t1',convertedfns);

    dcmhdrfn = fullfile(subjpath,'t1_dicom_header.mat');
    save(dcmhdrfn,'dcmhdr');
    aap = aas_desc_outputs(aap,i,'t1_dicom_header',dcmhdrfn);


    case 'checkrequirements'

    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end

