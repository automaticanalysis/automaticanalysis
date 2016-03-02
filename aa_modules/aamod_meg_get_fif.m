function [aap, resp] = aamod_meg_get_fif(aap,task,subj,sess)

resp='';

switch task
    case 'report'
        
    case 'doit'
        %% Locate series
        [d, megser] = aas_get_series(aap,'meg',subj,sess);
        
        %% Search file
        srcdir = meg_findvol(aap,aap.acq_details.subjects(subj).megname{d},'fp');
        megfile = fullfile(srcdir,megser);
        
        if ~exist(megfile,'file') % try as empty_room
            srcdir = meg_findvol(aap,megser,'fp');
            megfile = spm_select('FPList',srcdir,'.*fif');
            megfile = deblank(megfile(1,:)); % in case there are more, select the first
            aap.acq_details.subjects(subj).megseriesnumbers{d}{sess} = 'empty_room.fif';
        end
        
        if ~exist(megfile,'file')
            aas_log(aap,1,sprintf('ERROR: Subject %s has no session %s!',...
                aap.acq_details.subjects(subj).megname{d},...
                megser));
        end
        
        %% Copy file
        sessdir = aas_getsesspath(aap,subj,sess);
        copyfile(megfile,fullfile(sessdir,megser));
        
        %% Describe outputs
        aap=aas_desc_outputs(aap,subj,sess,'meg',...
            fullfile(sessdir,aap.acq_details.subjects(subj).megseriesnumbers{sess}));
end