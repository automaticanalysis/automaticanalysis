function [aap, resp] = aamod_meeg_maxfilt(aap,task,subj)

resp='';

switch task
    case 'report'
        
        aap = aas_report_add(aap,subj,'<table><tr>');
        
        for sess = aap.acq_details.selected_sessions
            aap = aas_report_add(aap,subj,'<td>');
            aap = aas_report_add(aap,subj,['<h3>Session: ' aas_getsessdesc(aap,subj,sess) '</h3>']);
            
            [sesspath, fstem] = fileparts(aas_getfiles_bystream(aap,'meeg_session',[subj,sess],'meeg','output'));
            logfname  = fullfile(sesspath,sprintf('%s.log',fstem));
            aap = aas_report_add(aap,subj,'<table><tr><td>');
            
            %% 1. Bad channels
            [junk, txt] = unix(sprintf('grep "bad channels (" %s',logfname));
            if ~isempty(txt)
                txt = textscan(txt,'%s','delimiter','\n'); txt = txt{1}{1};
                aap = aas_report_add(aap,subj,'<h4>Bad channels:</h4>');
                aap = aas_report_add(aap,subj,sprintf('<h5>%s</h5>',txt));
            end
            
            %% 2. Movement (courtesy to Jason Taylor and Rik Henson)
            [pth,fstem] = spm_fileparts(logfname);
            movfname = fullfile(pth,['diagnostic_aamod_meeg_maxfilt_' fstem 'movpar.txt']);
            figfname = fullfile(pth,['diagnostic_aamod_meeg_maxfilt_' fstem 'movpar.jpg']);
            
            % Parse log file for movement params:
            if ~exist(movfname,'file')
                unix(sprintf('sed -n -e ''/^#t/ p'' %s  | sed -e ''s/[^0-9,.]*//g'' -e ''s/#e = //g'' > %s',logfname,movfname));
            end
            
            % Load data:
            mov = [];
            d = dir(movfname);
            if d.bytes
                mov = dlmread(movfname,',');
            end
            
            % Check whether valid:
            aap = aas_report_add(aap,subj,'<h4>Movement:</h4>');
            if numel(mov)==0
                aap = aas_report_add(aap,subj,sprintf('<h5>WARNING: No movement data! - Was HPI estimation and movement compensation run in MF call that created %s?</h5>',movfname));
            else % plot
                if ~exist(figfname,'file')
                    try    fig = spm_figure('FindWin'); clf;
                    catch, fig = figure('color','w','paperpositionmode','auto');
                    end
                    time=mov(:,1);
                    subplot(2,1,1);
                    plot(time,mov(:,3),'g'); % goodness of fit
                    legend('          gof(0:1)','Location','NorthEastOutside');
                    title(basename(movfname),'Interpreter','none');
                    
                    subplot(2,1,2); hold on
                    plot(time,mov(:,2),'r'); % error
                    plot(time,mov(:,4),'k'); % velocity
                    plot(time,mov(:,5),'b'); % rotation
                    plot(time,mov(:,6),'c'); % translation
                    legend({'error(cm)','velocity(cm/s)','rotation(rad/s)','translation(cm)'},'Location','NorthEastOutside')
                    title(basename(movfname),'Interpreter','none');
                    
                    set(fig,'Renderer','zbuffer');
                    print(fig,'-djpeg',figfname);
                    aas_log(aap,false,sprintf('- Figure saved to %s\n',figfname));
                    close(fig)
                end
                
                % Return maximum translation, rotation:
                md = max(mov(:,6))-min(mov(:,6));
                mr = max(mov(:,5))-min(mov(:,5));
                
                aap = aas_report_addimage(aap,subj,figfname);
                aap = aas_report_add(aap,subj,sprintf('<h5>Max translation=%gcm. Max rotation=%grad/s.</h5>',md,mr));
            end
            
            %% 3. Transdef
            if ~isempty(aap.tasklist.currenttask.settings.transform)
                [sesspath, fstem] = fileparts(aas_getfiles_bystream(aap,'meeg_session',[subj,sess],'trans_meeg','output'));
                logtrfname  = fullfile(sesspath,sprintf('%s.log',fstem));
                
                aap = aas_report_add(aap,subj,'<h4>Transformation:</h4>');
                [junk, txt] = unix(sprintf('grep Rotation %s',logtrfname)); txt = txt(1:end-1);
                aap = aas_report_add(aap,subj,[fstem ': ' txt '<br>']);
                [junk, txt] = unix(sprintf('grep Position %s',logtrfname)); txt = txt(1:end-1);
                aap = aas_report_add(aap,subj,[fstem ': ' txt '<br>']);
            end
            
            aap = aas_report_add(aap,subj,'</td></tr></table>');
            
            aap = aas_report_add(aap,subj,'</td>');
        end
        
        aap = aas_report_add(aap,subj,'</tr></table>');
    case 'doit'
        
        %% If session for HPI is specified, bring forward
        sessions = aap.acq_details.selected_sessions;
        HPIsess = aas_getsetting(aap,'hpi.session');
        if ~isempty(HPIsess)
            HPIsess = cell_index({aap.acq_details.meeg_sessions(aap.acq_details.selected_sessions).name},HPIsess);
            sessions = [HPIsess sessions];
            [junk, ind] = unique(sessions);
            sessions = sessions(sort(ind));
        end
        
        for sess = sessions
                %% Initialise
                sesspath = aas_getsesspath(aap,subj,sess);
                
                % clear previous diagnostics
                delete(fullfile(sesspath,'diagnostic_aamod_meeg_maxfilt_*'));
                
                instream = aas_getstreams(aap,'input'); instream = instream{1};
                infname = aas_getfiles_bystream(aap,'meeg_session',[subj sess],instream);
                outfname = fullfile(sesspath,['mf2pt2_' basename(infname) '.fif']); % specifying output filestem
                delete(fullfile(sesspath,'*mf2pt2_*'));
                [pth, fstem, ext] = fileparts(outfname);
                
                isEmptyRoom = strcmp(aap.acq_details.meeg_sessions(sess).name,'empty_room');
                
                %% Sphere fit
                spherefit = [];
                if ~isEmptyRoom % if empty_room session, skip ...
                    addpath(fullfile(aap.directory_conventions.neuromagdir,'meg_pd_1.2'));
                    try spherefit = meg_fit_sphere_rik(aap,infname,outfname);
                    catch E
                        aas_log(aap,1,sprintf('ERROR: Spherefit failed for %s: %s!\nERROR: No frame/origin will be applied!',aas_getsessdesc(aap,subj,sess),E.message))
                    end;
                    rmpath(fullfile(aap.directory_conventions.neuromagdir,'meg_pd_1.2'));
                else
                    if ~isempty(HPIsess) % ... and load from HPIsess
                        meegfname = aas_getfiles_bystream(aap,'meeg_session',[subj,HPIsess],'meeg','output');
                        spherefit = load(spm_file(meegfname,'suffix','_sphere_fit','ext','txt'),'-ASCII','spherefit');
                    else
                        aas_log(aap,1,sprintf('WARNING: Session for HPI not defined!\nWARNING: HPI is not perfomed for %s!',aas_getsessdesc(aap,subj,sess)))
                    end
                end
                
                %% Maxfilt
                orgcmd = '';
                hpicmd = '';
                stcmd  = '';
                trcmd_par = '';
                
                % Origin and frame
                if ~isempty(spherefit)
                    orgcmd = sprintf(' -frame head -origin %g %g %g',spherefit(1),spherefit(2),spherefit(3)');
                end
                
                % Autobad throughout (assumes 15mins=900s)
                badstr  = sprintf(' -autobad %d -badlimit %d',...
                    aas_getsetting(aap,'autobad.interval',sess),...
                    aas_getsetting(aap,'autobad.badlimit',sess));
                
                % SSS with ST:
                if aap.tasklist.currenttask.settings.sss.run
                    stcmd    = sprintf(' -st %d -corr %g ',...
                        aas_getsetting(aap,'sss.window',sess),...
                        aas_getsetting(aap,'sss.corr',sess));
                    
                else
                    stcmd = ''; % To avoid jumps between end of 10s buffer, as warned in manual
                end
                
                % Do movecomp?
                if strcmp(aap.tasklist.currenttask.settings.sss.mvcomp, 'on') || isempty(aap.tasklist.currenttask.settings.sss.mvcomp)
                    mvcmd = ' -movecomp inter'; % areguments added to the maxfilter unix command
                    mvcomp_out = char(fullfile(sesspath,[fstem '_headpoints.txt']),...
                        fullfile(sesspath,[fstem '_headposition.pos'])); % used later when checking md5 checksums
                elseif strcmp(aap.tasklist.currenttask.settings.sss.mvcomp, 'off')
                    mvcmd = '';
                    mvcomp_out = char(fullfile(sesspath,[fstem '_headpoints.txt'])); % used later when checking md5 checksums (no headpos in this case)
                else
                    aas_log(aap,true,sprintf('ERROR: invalid option for sss.movecomp valid inputs = on|off'));
                end
                
                % HPI estimation and movement compensation
                if ~isEmptyRoom
                    hpicmd = sprintf(' -linefreq 50 -hpistep %d -hpisubt %s -hpicons -hp %s',...
                        aas_getsetting(aap,'hpi.step',sess),...
                        aas_getsetting(aap,'hpi.subt',sess),...
                        fullfile(sesspath,[fstem '_headposition.pos']));
                else
                    hpicmd = ' -linefreq 50';
                end
                
                % Preparing names of output and log files
                outfname  = fullfile(sesspath,sprintf('%s.fif',fstem));
                logfname  = fullfile(sesspath,sprintf('%s.log',fstem));
                
                skipstr = '';
                
                mfcall = fullfile(aap.directory_conventions.neuromagdir,'bin','util','maxfilter-2.2.12');
                
                % calibration files
                fnamectc = fullfile(aap.directory_conventions.neuromagdir,'databases','ctc','ct_sparse.fif');
                fnamecal = fullfile(aap.directory_conventions.neuromagdir,'databases','sss','sss_cal.dat');
                calibdir = aas_getsetting(aap,'calibrationdirectory');
                if ~isempty(calibdir)
                    if ~exist(calibdir,'dir')
                        aas_log(aap,true,sprintf('ERROR: Calibration directory %s not exist!',calibdir));
                    else
                        if exist(fullfile(calibdir,'ct_sparse.fif'),'file'), fnamectc = fullfile(calibdir,'ct_sparse.fif'); end
                        if exist(fullfile(calibdir,'sss_cal.dat'),'file'), fnamecal = fullfile(calibdir,'sss_cal.dat'); end
                    end
                end
                
                % Assembling MF command
                mfcmd_rest=[
                    mfcall ' -f ' infname ' -o ' outfname,...
                    ' -ctc ' fnamectc,...
                    ' -cal ' fnamecal,...
                    ' ', skipstr, badstr, orgcmd, stcmd, hpicmd, mvcmd, trcmd_par ' -force -v | tee ' logfname
                    ];
                aas_log(aap,false,mfcmd_rest);
                
                % Executing MF
                [status, maxres] = unix(mfcmd_rest); % this stops screen-dumping?
                if status
                    aas_log(aap,true,'ERROR: MaxFilter failed!')
                else
                    aas_log(aap,false,maxres);
                end
                
                %% Trans (so that origin not same as SSS expansion origin above)
                if ~isempty(spherefit) && ~isempty(aap.tasklist.currenttask.settings.transform)
                    outtrfname = outfname;
                    outtrpfx = '';
                    for t = 1:numel(aap.tasklist.currenttask.settings.transform)
                        if iscell(aap.tasklist.currenttask.settings.transform)
                            trans = aap.tasklist.currenttask.settings.transform{t};
                        else
                            trans = aap.tasklist.currenttask.settings.transform(t);
                        end
                        if ischar(trans) % fif file
                            ref_str = trans;
                            outtrpfx    = ['trans' basename(trans) '_' outtrpfx];
                        elseif isnumeric(trans)
                            if trans
                                if trans < 0 % selected_session number
                                    trans = aap.acq_details.selected_sessions(-trans);
                                end
                                if trans == sess, continue; end % do not trans to itself
                                if isEmptyRoom && (trans == HPIsess), continue; end % do not trans empty_room to HPIsess
                                ref_str = aas_getfiles_bystream(aas_setcurrenttask(aap,cell_index({aap.tasklist.main.module.name},'aamod_meeg_get_fif')),... % raw data
                                    'meeg_session',[subj,trans],instream);
                                trcmd_par = sprintf(' -trans %s ',ref_str);
                                outtrpfx    = ['trans' aap.acq_details.meeg_sessions(trans).name '_' outtrpfx];
                            else % 0
                                trcmd_par = sprintf(' -trans default -origin %g %g %g -frame head ',spherefit(1),spherefit(2)-13,spherefit(3)+6);
                                outtrpfx    = ['transdef_'  outtrpfx];
                            end
                        else
                            aas_log(aap,1,'Trans reference: Unrecognised option!');
                        end
                        intrfname   = outtrfname;
                        outtrfname  = fullfile(sesspath,sprintf('%s%s.fif',outtrpfx,fstem));
                        logtrfname  = fullfile(sesspath,sprintf('%s%s.log',outtrpfx,fstem));
                        
                        % Assembling MF command
                        mfcmd_rest=[
                            mfcall ' -f ' intrfname ' -o ' outtrfname,...
                            trcmd_par ' -force -v | tee ' logtrfname
                            ];
                        aas_log(aap,false,mfcmd_rest);
                        
                        % Executing MF
                        [status, maxres] = unix(mfcmd_rest); % this stops screen-dumping?
                        if status ~= 0
                            aas_log(aap,true,'ERROR: Trans MaxFilter failed!')
                        else
                            aas_log(aap,false,maxres);
                        end
                    end
                end
                
                outs{1} = outfname; logs{1} = logfname;
                if ~isempty(spherefit) && ~isempty(aap.tasklist.currenttask.settings.transform)
                    outs{2} = outtrfname; logs{2} = logtrfname;
                end
                
                %% Downsampling
                if ~isempty(aas_getsetting(aap,'downsampling'))
                    for i = 1:numel(outs)
                        outdsfname = spm_file(outs{i},'prefix','ds_');
                        mfcmd_ds=[
                            mfcall ' -f ' outs{i} ' -o ' outdsfname,...
                            ' -ds ' num2str(aas_getsetting(aap,'downsampling')) ' -force -v | tee --append ' logs{i}
                            ];
                        aas_log(aap,false,mfcmd_ds);
                        
                        % Executing MF
                        [status, maxres] = unix(mfcmd_ds); % this stops screen-dumping?
                        if status ~= 0
                            aas_log(aap,true,'ERROR: Trans MaxFilter failed!')
                        else
                            aas_log(aap,false,maxres);
                        end
                        outs{i} = outdsfname;
                    end
                end
                
                %% Outputs
                outstream = aas_getstreams(aap,'output');
                [junk,ind] = min(cellfun(@(x) numel(x), outstream)) % shortest
                
                aap=aas_desc_outputs(aap,subj,sess,outstream{ind},outs{1});
                
                if ~isEmptyRoom
                    aap=aas_desc_outputs(aap,subj,sess,outstream{contains(outstream,'head')},mvcomp_out);
                end
                
                if ~isempty(spherefit) && ~isempty(aap.tasklist.currenttask.settings.transform)
                    aap=aas_desc_outputs(aap,subj,sess,outstream{contains(outstream,'trans')},outs{2});
                end
        end
end
