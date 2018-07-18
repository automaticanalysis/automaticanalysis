% AA module - Converts fieldmaps maps to NIFTI format
% Rhodri Cusack MRC CBU Cambridge Nov 2005

function [aap,resp]=aamod_convert_fieldmaps(aap,task,subj,sess)

resp='';

switch task
    case 'report'
    case 'doit'
        domain = aap.tasklist.currenttask.domain;
        
        [aap, convertedfns, dcmhdr]=aas_convertseries_fromstream(aap,domain,[subj,sess],'dicom_fieldmap');
        
        sesspath=aas_getsesspath(aap,subj,sess);
        
        % Restructure outputs
        outstream = char(vertcat(convertedfns{:}));
        aap=aas_desc_outputs(aap,domain,[subj,sess],'fieldmap',outstream);

        if iscell(dcmhdr{1}) % multiple series --> selects the one with the most scans (it should contain all TEs)
            [junk,ind] = max(cellfun(@(x) numel(x), dcmhdr));
            dcmhdr = dcmhdr{ind};
        end
        dcmhdrfn=fullfile(sesspath,'fieldmap_dicom_header.mat');
        save(dcmhdrfn,'dcmhdr');
        aap=aas_desc_outputs(aap,domain,[subj,sess],'fieldmap_dicom_header',dcmhdrfn);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,true,sprintf('Unknown task %s',task));
end