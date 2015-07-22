function [aap, resp] = aamod_meg_get_fif(aap,task,subj,sess)

resp='';

switch task
    case 'report'
        
    case 'doit'
        %% Search file
        srcdir = meg_findvol(aap,aap.acq_details.subjects(subj).megname,'fp');
        megfile = fullfile(srcdir,aap.acq_details.subjects(subj).megseriesnumbers{sess});
        
        if ~exist(megfile,'file') % try as empty_room
            srcdir = meg_findvol(aap,aap.acq_details.subjects(subj).megseriesnumbers{sess},'fp');
            megfile = spm_select('FPList',srcdir,'.*fif');
            aap.acq_details.subjects(subj).megseriesnumbers{sess} = 'empty_room.fif';
        end
        
        if ~exist(megfile,'file')
            aas_log(aap,1,sprintf('ERROR: Subject %s has no session %s!',...
                aap.acq_details.subjects(subj).megname,...
                aap.acq_details.subjects(subj).megseriesnumbers{sess}));
        end
        
        %% Copy file
        sessdir = aas_getsesspath(aap,subj,sess);
        copyfile(megfile,fullfile(sessdir,aap.acq_details.subjects(subj).megseriesnumbers{sess}));
        
        %% Describe outputs
        aap=aas_desc_outputs(aap,subj,sess,'meg',...
            fullfile(sessdir,aap.acq_details.subjects(subj).megseriesnumbers{sess}));
end