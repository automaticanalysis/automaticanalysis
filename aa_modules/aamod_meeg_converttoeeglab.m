function [aap, resp] = aamod_meeg_converttoeeglab(aap,task,subj,sess)

resp='';

switch task
    case 'report'
        for fn = cellstr(spm_select('FPList',aas_getsesspath(aap,subj,sess),'^diagnostic_.*jpg$'))'
            aap = aas_report_add(aap,subj,'<table><tr><td>');
            aap=aas_report_addimage(aap,subj,fn{1});
            aap = aas_report_add(aap,subj,'</td></tr></table>');
        end
    case 'doit'
        infname = cellstr(aas_getfiles_bystream(aap,'meeg_session',[subj sess],'meeg'));
        
        [junk, EL] = aas_cache_get(aap,'eeglab');
        EL.load;
        
        % read data
        EEG = [];
        for i = 1:numel(infname)
            try
                EEG = pop_fileio(infname{i});
                [res,o1,o2] = evalc('eeg_checkset(EEG)');
                if isempty(res), break;
                else, throw(res); end
            catch err
                E(i) = err;
            end
        end
        if isempty(EEG)
            for i = 1:numel(E)
                aas_log(aap,i==numel(E),sprintf('ERROR: reading %s - %s',infname{i}, E(i).message));
            end
        end
        
        % channel layout
        EEG = pop_chanedit(EEG,'lookup',aas_getfiles_bystream(aap,'channellayout'));
        
        % remove channel
        if ~isempty(aas_getsetting(aap,'removechannel'))
            chns = strsplit(aas_getsetting(aap,'removechannel'),':');
            EEG = pop_select(EEG, 'nochannel', cellfun(@(x) find(strcmp({EEG.chanlocs.labels}, x)), chns));
        end
        
        % downsample
        if ~isempty(aas_getsetting(aap,'downsample'))
            sRate = aas_getsetting(aap,'downsample');
            if sRate ~= EEG.srate, EEG = pop_resample( EEG, aas_getsetting(aap,'downsample')); end
        end
        
        % edit
        % - specify operations
        toEditsetting = aas_getsetting(aap,'toEdit');
        toEditsubj = toEditsetting(...
            cellfun(@(x) any(strcmp(x,aas_getsubjname(aap,subj))),{toEditsetting.subject}) | ...
            strcmp({toEditsetting.subject},'*')...
            );        
        toEdit = struct('type',{},'operation',{});
        for s = 1:numel(toEditsubj)
            sessnames = regexp(toEditsubj(s).session,':','split');
            if any(strcmp(sessnames,aas_getsessname(aap,sess))) || sessnames{1} == '*'
                toEdit = horzcat(toEdit,toEditsubj(s).event);
            end
        end
        
        % - do it
        if ~isempty(toEdit)
            for e = toEdit
                if ischar(e.type)
                    ind = ~cellfun(@isempty, regexp({EEG.event.type},e.type));
                elseif isnumeric(e.type)
                    ind = e.type;
                end
                op = strsplit(e.operation,':');
                if ~any(ind) && ~strcmp(op{1},'insert'), continue; end
                switch op{1}
                    case 'remove'
                        EEG.event(ind) = [];
                        EEG.urevent(ind) = [];
                    case 'keep'
                        EEG.event = EEG.event(ind);
                        EEG.urevent = EEG.urevent(ind);
                    case 'rename'
                        for i = find(ind)
                            EEG.event(i).type = op{2};
                            EEG.urevent(i).type = op{2};
                        end
                    case 'unique'
                        ex = [];
                        switch op{2}
                            case 'first'
                                for i = 2:numel(ind)
                                    if ind(i) && ind(i-1), ex(end+1) = i; end
                                end
                            case 'last'
                                for i = 1:numel(ind)-1
                                    if ind(i) && ind(i+1), ex(end+1) = i; end
                                end
                        end
                        EEG.event(ex) = [];
                        EEG.urevent(ex) = [];
                    case 'iterate'
                        ind = cumsum(ind).*ind;
                        for i = find(ind)
                            EEG.event(i).type = sprintf('%s%02d',EEG.event(i).type,ind(i));
                            EEG.urevent(i).type = sprintf('%s%02d',EEG.urevent(i).type,ind(i));
                        end
                    case 'insert'
                        loc = str2num(op{2});
                        newE = EEG.event(loc);
                        for i = 1:numel(newE)
                            newE(i).type = e.type;
                        end
                        events = EEG.event(1:loc(1)-1);
                        for i = 1:numel(loc)-1
                            events = [events newE(i) EEG.event(loc(i):loc(i+1)-1)];
                        end
                        if isempty(i), i = 0; end
                        events = [events newE(i+1) EEG.event(loc(i+1):end)];
                        EEG.event = events;
                        EEG.urevent = rmfield(events,'urevent');
                    case 'ignorebefore'
                        EEG = pop_select(EEG,'nopoint',[0 EEG.event(ind(1)).latency-1]);
                        beInd = find(strcmp({EEG.event.type},'boundary'),1,'first');
                        samplecorr = EEG.event(beInd).duration;
                        EEG.event(beInd) = [];
                        ureindcorr = EEG.event(1).urevent -1;
                        
                        % adjust events                        
                        for i = 1:numel(EEG.event)
                            EEG.event(i).urevent = EEG.event(i).urevent - ureindcorr;
                        end
                        EEG.urevent(1:ureindcorr) = [];
                        
                        % adjust time
                        for i = 1:numel(EEG.urevent)
                            EEG.urevent(i).latency = EEG.urevent(i).latency - samplecorr;
                        end
                    case 'ignoreafter'
                        EEG = pop_select(EEG,'nopoint',[EEG.event(ind(end)).latency EEG.pnts]);
                        ureindcorr = EEG.event(end).urevent;
                        
                        % adjust events
                        EEG.urevent(ureindcorr+1:end) = [];
                    otherwise
                        aas_log(aap,false,sprintf('Operation %s not yet implemented',op{1}));
                end
                % update events
                for i = 1:numel(EEG.event)
                    urind = find(strcmp({EEG.urevent.type},EEG.event(i).type) & [EEG.urevent.latency]==EEG.event(i).latency);
                    EEG.event(i).urevent = urind(1);
                end
            end
        end

        % diagnostics
        diagpath = fullfile(aas_getsesspath(aap,subj,sess),['diagnostic_' mfilename '_raw.jpg']);
        meeg_diagnostics_continuous(EEG,aas_getsetting(aap,'diagnostics'),'Raw',diagpath);
                
        fnameroot = sprintf('eeg_%s',aas_getsubjname(aap,subj));
        while ~isempty(spm_select('List',aas_getsesspath(aap,subj,sess),[fnameroot '.*']))
            fnameroot = ['aa' fnameroot];
        end
        
        pop_saveset(EEG,'filepath',aas_getsesspath(aap,subj,sess),'filename',fnameroot);
        outfname = spm_select('FPList',aas_getsesspath(aap,subj,sess),[fnameroot '.*']);
        
        EL.unload;
        
        %% Describe outputs
        aap = aas_desc_outputs(aap,subj,sess,'meeg',outfname);
    case 'checkrequirements'
        if ~aas_cache_get(aap,'eeglab'), aas_log(aap,true,'EEGLAB is not found'); end
end