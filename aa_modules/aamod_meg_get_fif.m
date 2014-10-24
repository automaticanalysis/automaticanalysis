function [aap, resp] = aamod_meg_get_fif(aap,task,subj,sess)

resp='';

switch task
    case 'report'
        
    case 'doit'
        sessdir = aas_getsesspath(aap,subj,sess);
        tmpaap = aap;
        tmpaap.directory_conventions.megsubjectoutputformat = '%s';
        srcdir = meg_findvol(tmpaap,tmpaap.acq_details.subjects(subj).megname,1);
        megfile = fullfile(srcdir,aap.acq_details.subjects(subj).megseriesnumbers{sess});
        copyfile(megfile,fullfile(sessdir,aap.acq_details.subjects(subj).megseriesnumbers{sess}));
        aap=aas_desc_outputs(aap,subj,sess,'meg',...
            fullfile(sessdir,aap.acq_details.subjects(subj).megseriesnumbers{sess}));
end