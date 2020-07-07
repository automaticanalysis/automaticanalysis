function [aap, resp] = aamod_meeg_get_file(aap,task,subj,sess)

resp='';

switch task
    case 'report'
        
    case 'doit'
        %% Locate series
        [d, meegser] = aas_get_series(aap,'MEEG',subj,sess);
        
        %% Search file
        if exist(meegser,'file') % full path
            srcdir = spm_file(meegser,'path');
            meegser = spm_file(meegser,'filename');
        else % relative to subject's folder
            srcdir = meeg_findvol(aap,aap.acq_details.subjects(subj).meegname{d},'fullpath',true); % subject's folder
            srcdir = fullfile(srcdir,spm_file(meegser,'path')); % any subfolder
            meegser = spm_file(meegser,'filename');
        end
                
        if ~exist(fullfile(srcdir,meegser),'file') % try as empty_room
            srcdir = meeg_findvol(aap,meegser,'fullpath',true);
            meegser = spm_select('List',srcdir,'.*fif');
            meegser = deblank(meegser(1,:)); % in case there are more, select the first
            aap.acq_details.subjects(subj).meegseriesnumbers{d}{sess} = 'empty_room.fif';
        end
        
        meegfile = fullfile(srcdir,meegser);
        
        if ~exist(meegfile,'file')
            aas_log(aap,1,sprintf('ERROR: Subject %s has no session %s!',...
                aap.acq_details.subjects(subj).meegname{d},...
                meegser));
        end
        
        %% Copy file
        sessdir = aas_getsesspath(aap,subj,sess);
        
        switch spm_file(meegser,'ext')
            case 'fif' % MEG Neuromag
                copyfile(meegfile,fullfile(sessdir,meegser),'f'); % DP changed to force file overwrite
            case 'vhdr' % EEG BrainVision
                baseFn = spm_file(meegser,'basename');
                meegser = {};
                for f = cellstr(spm_select('FPList',srcdir,['^' baseFn '.*']))'
                    meegser{end+1} = spm_file(f{1},'filename');
                    copyfile(f{1},fullfile(sessdir,meegser{end}),'f'); % DP changed to force file overwrite
                end
        end

        %% Describe outputs
        aap=aas_desc_outputs(aap,subj,sess,'meeg',...
            fullfile(sessdir,meegser));
end