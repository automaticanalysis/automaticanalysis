function [aap, resp] = aamod_meg_maxfilt(aap,task,subj,sess)

resp='';

switch task
    case 'report'
        
    case 'doit'
        %% Initialise
        sessdir = aas_getsesspath(aap,subj,sess);
        infname = aas_getfiles_bystream(aap,'meg_session',[subj sess],'meg');        
        outfname = fullfile(sessdir,['mf2pt2_' basename(infname) '.fif']); % specifying output filestem
        delete(fullfile(sessdir,'*mf2pt2_*'));
        [pth, fstem, ext] = fileparts(outfname);
        
        %% Sphere fit
        addpath(fullfile(aap.directory_conventions.neuromagdir,'meg_pd_1.2'));
        spherefit = meg_fit_sphere_rik(aap,infname,outfname);
        rmpath(fullfile(aap.directory_conventions.neuromagdir,'meg_pd_1.2'));

        %% Maxfilt
        orgcmd = '';
        hpicmd = '';
        stcmd  = '';
        trcmd_par = '';
        
        % Origin and frame
        orgcmd = sprintf(' -frame head -origin %g %g %g',spherefit(1),spherefit(2),spherefit(3)');

        % Autobad throughout (assumes 15mins=900s)
        badstr  = sprintf(' -autobad %d -badlimit %d',...
            aap.tasklist.currenttask.settings.autobad.interval,...
            aap.tasklist.currenttask.settings.autobad.badlimit);
        
        % SSS with ST:
        if aap.tasklist.currenttask.settings.sss.run
            stcmd    = sprintf(' -st %d -corr %g ',...
                aap.tasklist.currenttask.settings.sss.window,...
                aap.tasklist.currenttask.settings.sss.corr);
        else
            stcmd = ''; % To avoid jumps between end of 10s buffer, as warned in manual
        end

        % HPI estimation and movement compensation
        hpicmd = sprintf(' -linefreq 50 -hpistep %d -hpisubt %s -hpicons -movecomp inter -hp %s',...
            aap.tasklist.currenttask.settings.hpi.step,...
            aap.tasklist.currenttask.settings.hpi.subt,...
            fullfile(sessdir,[fstem '_headposition.pos']));

        % Preparing names of output and log files
        outpfx    = '';
        outfname  = fullfile(sessdir,sprintf('%s%s.fif',outpfx,fstem));
        logfname  = fullfile(sessdir,sprintf('%s%s.log',outpfx,fstem));
        
        skipstr = '';
        
        mfcall = fullfile(aap.directory_conventions.neuromagdir,'bin','util','maxfilter-2.2.12');
        
        % Assembling MF command
        mfcmd_rest=[
            mfcall ' -f ' infname ' -o ' outfname,...
            ['	 -ctc ' fullfile(aap.directory_conventions.neuromagdir,'databases','ctc','ct_sparse.fif')] ' ',...
            ['	 -cal ' fullfile(aap.directory_conventions.neuromagdir,'databases','sss','sss_cal.dat')] ' ',...
            skipstr, badstr, orgcmd, stcmd, hpicmd, trcmd_par ' -v | tee ' logfname
            ];
        disp(mfcmd_rest);

        % Executing MF
        [status, maxres] = unix(mfcmd_rest); % this stops screen-dumping?
        if status
            aas_log(aap,1,'MaxFilter failed!')
        else
            disp(maxres);
        end
        
        %% Trans (so that origin not same as SSS expansion origin above)
        if ~isempty(aap.tasklist.currenttask.settings.transform)
            if aap.tasklist.currenttask.settings.transform % session number
                ref_str = aas_getfiles_bystream(aas_setcurrenttask(aap,1),... % raw data
                    'meg_session',[subj,aap.tasklist.currenttask.settings.transform],'meg');
                outpfx    = ['trans' aap.acq_details.meg_sessions(aap.tasklist.currenttask.settings.transform),name '_'];
            else % 0
                ref_str = 'default' ;
                outpfx    = 'transdef_';
            end
            trcmd_par = sprintf(' -trans %s -origin %g %g %g -frame head ',ref_str,spherefit(1),spherefit(2)-13,spherefit(3)+6);
            infname   = outfname;
            outtrfname  = fullfile(sessdir,sprintf('%s%s.fif',outpfx,fstem));
            logtrfname  = fullfile(sessdir,sprintf('%s%s.log',outpfx,fstem));

            % Assembling MF command
            mfcmd_rest=[
                mfcall ' -f ' infname ' -o ' outtrfname,...
                trcmd_par ' -force -v | tee ' logtrfname
                ];
            disp(mfcmd_rest);

            % Executing MF
            [status, maxres] = unix(mfcmd_rest); % this stops screen-dumping?
            if status ~= 0
                aas_log(aap,1,'Trans MaxFilter failed!')
            else
                disp(maxres);
            end
        end
        
        %% Outputs
        aap=aas_desc_outputs(aap,subj,sess,'meg_head',...
            char(fullfile(sessdir,[fstem '_headpoints.txt']),...
            fullfile(sessdir,[fstem '_headposition.pos'])));        
        aap=aas_desc_outputs(aap,subj,sess,'meg',outfname);        
        if aap.tasklist.currenttask.settings.transform
            aap=aas_desc_outputs(aap,subj,sess,'norm_meg',outtrfname);
        end
end