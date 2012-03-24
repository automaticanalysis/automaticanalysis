function [aap resp]=aamod_emeg_inversion(varargin)
% Do inversions. Now invert each condition and time window seperately to
% avoid signal leakage.
% Danny Mitchell 02/04/08; 07/11/08

%% check task settings, subject, block etc
[aap subblock doit resp settings]=aa_emeg_checktasksettings(mfilename('fullpath'),varargin);

%% CHECK REQUIREMENTS: load spreadsheet containing epoching instructions
epochfile=spm_select('FPList',fullfile(aas_getsesspath(aap,subblock{1},subblock{2}),'events'),'(E|e)pochs.*.xls');
if exist(epochfile,'file')~=2; epochfile=spm_select('FPList',fullfile(aas_getsubjpath(aap,subblock{1}),'events'),'(E|e)pochs.*\.xls'); end
if exist(epochfile,'file')~=2; epochfile=spm_select('FPList',fullfile(aap.acq_details.root,'events'),'(E|e)pochs.*.xls'); end
if exist(epochfile,'file')~=2; error('aa:EpochFileError', '\n Failed to find Epochs.*.xls, required for epoching continuous data. \n This should be placed in "events" folder at either study, subject, or session levels, with lower levels overriding higher ones. \n See examples/Epochs.xls for an example of the required format. \n'); end

if ~doit; return; end
try settings=rmfield(settings,'specialrequirements'); catch; end

%% find files and decide whether to run task;
files=aas_emeg_findfiles(aap,settings.InputFilter,subblock);
if isempty(files); aas_log(aap,1,sprintf('\nFound no data! (Input filter is %s)\n',settings.InputFilter)); end
clear varargin

%% load event specifications from spreadsheet
warning off all; [Numeric,Txt]=xlsread(epochfile); warning on all

%% for each file excluding EEG or CDA
for f=1:length(files);

    %% skip EEG or CDA files
    if ~isempty(regexp(files{f},'-eeg','ONCE')); continue; end
    if ~isempty(regexp(files{f},'_CDA_','ONCE')); continue; end

    %% load MEG header
    [pth hdrfn]=fileparts(files{f});
    hdrfn=fullfile(pth,[hdrfn '.mat']);
    try load(hdrfn);
    catch
        fprintf('\nFailed to load file %s.\nWill try again in 60s in case there is an access conflict....\n',hdrfn);
        pause(60)
        try
            load(hdrfn);
        catch
            hdrfn=strrep(hdrfn,'.mat','.bck');
            if exist(hdrfn,'file')
                fprintf('Failed again. Trying to load backup...\n')
                try load(hdrfn,'-MAT')
                catch aas_log(aap,1,sprintf('Failed again. File %s may be corrupt?',hdrfn(1:end-4)));
                end
            end
            %            fprintf('\nDebug %s\n',mfilename); keyboard
        end
    end

    %% delete unnecessary fields to save memory
    redir=sprintf('See %s',D.fname(2:end));
    if isfield(D,'thresholds') && ~strcmp(D.fname(1),'a'); D.thresholds=redir; end
    if isfield(D,'filter') && ~strcmp(D.fname(1),'f'); D.filter=redir; end
    if isfield(D,'ica') && ~strcmp(D.fname(1),'p'); D.ica=redir; end

    %% make directory for inversions
    if ~exist(fullfile(D.path,'inversions'),'dir')
        mkdir(fullfile(D.path,'inversions'))
    end

    %% load relevant events/contrasts
    try
        anchor=~cellfun('isempty',regexp(D.fname,Txt(1,:)));
        col=~cellfun('isempty',regexp(regexprep(Txt(2,anchor),{'[',']'},{'',''}),regexprep(Txt(2,:),{'[',']'},{'',''})));
    catch continue
    end

    %% get windows for time/frequency contrasts
    w=Txt(4:5,col);
    TimeWindows=str2num(w{1}); % matrix
    FrequencyWindows=str2num(w{2}); % matrix

    %% find columns containing events and zeroed contrasts
    n=Numeric(:,col);
    n(all(isnan(n')),:)=[]; % remove rows of NaNs; only happened for one subject, possibly due to formatting of .xls
    Ecol=find(size(n,1)-sum(n==0,1)==1); % relevant columns containing events
    Ccol=find(round(sum(n))==0); % relevant columns containing zeroed contrasts
    enames=Txt(1,col);

    %% determine sensor type(s)
    stype=regexp(files{f},'(-\w\w\w\w)\.','tokens');
    if isempty(stype); stype=''; % both sensor types
    else stype=stype{1}{1};
    end

    nonfuse=true; fuse=false;
    if isempty(stype); % data contains both mags and grads

        %% if data contains both mags and grads, decide whether to try
        %% fusion inversion and/or standard inversion
        % 1=fusion; 2=normal; 3=both
        if settings.ReweightSensorTypes==1; nonfuse=false; end
        if settings.ReweightSensorTypes~=2 % try fusion
            [pth nam ext]=fileparts(files{f});
            magfile=fullfile(pth,['s' nam '-mags' ext]);
            grdfile=fullfile(pth,['s' nam '-grds' ext]);
            if exist(magfile,'file') && exist(grdfile,'file')
                Dm=spm_eeg_ldata(magfile);
                Dg=spm_eeg_ldata(grdfile);
                try
                    Dm.inv{1}.forward.method;
                    Dg.inv{1}.forward.method;
                    fuse=true;
                catch aas_log(aap,true,sprintf('\nFound seperate files for mags & grds, but not forward models - NOT fusing.'))
                end
            else
                aas_log(aap,true,sprintf('\nFailed to find seperate file for mags & grds - NOT fusing.'))
            end;
        end
    end
    if ~(fuse||nonfuse); continue; end

    fprintf('\nProcessing %s\n',files{f});
    clear anchor doit epochfile hdrfn pth redir w

    %% find all forward model types in this file
    forwardmodels={};
    try
        for inv=1:length(D.inv)
            forwardmodels=[forwardmodels D.inv{inv}.forward.method];
        end
    catch aas_log(aap,false,'Check inversion');debugnow
    end
    forwardmodels=unique(forwardmodels);

    try D=rmfield(D,'zucosl');catch; end % relic: zeroed unsigned contrast of signed localisation
    % Seems invalid as it will pick up noise.
    % E.g. IPS equally active for conditions A and B.
    % Contrast will have expectation zero, but variance (varA + varB).
    % Absolute value of contrast probably greater for region activiated in
    % both conditions than region activated in neither.
    % Same applies to power of localisation
    try D=rmfield(D,'zscoul'); catch; end % relic: no longer zeroed
    try D=rmfield(D,'coul'); catch; end % relic
    try D=rmfield(D,'cosl'); catch; end % relic
    if ~isfield(D,'coep') || ~iscell(D.coep) || settings.Overwrite; D.coep={}; end % will be filled for each inversion type
    if ~isfield(D,'coeJ') || ~iscell(D.coeJ) || settings.Overwrite; D.coeJ={}; end % will be filled for each inversion type
    if ~isfield(D,'coip') || ~iscell(D.coip) || settings.Overwrite; D.coip={}; end
    if ~isfield(D,'cotp') || ~iscell(D.cotp) || settings.Overwrite; D.cotp={}; end

    %% do standard inversion then fusion inversion as requested
    for dofuse=0:1
        if ~dofuse && ~nonfuse; continue; end
        if dofuse && ~fuse; continue; end
        %% do this inversion for all types of forward model
        for fm=1:length(forwardmodels)
            for tw=1:size(TimeWindows,1)
                condinvs=zeros(size(D.events.types));
                for condnum=1:length(D.events.types)
                    cond=D.events.types(condnum);
                    %% Prepare inversion parameters
                    %     inverse.trials - D.events.types to invert
                    %     inverse.smooth - smoothness of source priors (0 to 1;
                    %     default 0.6 which, for a 4004 vertex mesh,
                    %     propogates spatial dependencies over 3-4 vertices that are
                    %     on average 6mm apart. See Friston et al. submitted 2007)
                    %     inverse.Np     - number of sparse priors per hemisphere (def. 256)
                    %     inverse.Nm     - maximum number of channel modes
                    %     inverse.NmTOL  - (exponentiated) threshold for spatial SVD (if Nm=[])
                    %     inverse.Nr     - maximum number of temporal modes
                    %     inverse.NrTOL  - (exponentiated) threshold for temporal SVD (if Nr=[])
                    %     inverse.type   - 'GS' Greedy search on MSPs
                    %                      'ARD' ARD search on MSPs
                    %                      'MSP' GS and ARD multiple sparse priors
                    %                      'LOR' LORETA-like model
                    %                      'IID' LORETA and minimum norm
                    %     inverse.xyz    - (n x 3) locations of spherical VOIs (def. [0 0 0])
                    %     inverse.rad    - radius (mm) of VOIs (def. 128)
                    %     inverse.lpf    - band-pass filter - low  frequency cutoff (Hz) (def. 1)
                    %     inverse.hpf    - band-pass filter - high frequency cutoff (Hz) (def. 256, but max 48 for MEG!)
                    %     inverse.Lap    - switch for Laplace transform
                    %     inverse.sdv    - standard devations of Gaussian temporal correlation (def. 4ms)
                    %     inverse.Han    - switch for Hanning window (def. 1)
                    %     inverse.Na     - number of most energetic dipoles (def. 1024)
                    %     inverse.woi    - time window for inversion ([start stop] in ms; def. all; NB not the window used for conditional contrasts later...)
                    % Evaluates:
                    %     inverse.M      - MAP projector (reduced)
                    %     inverse.J      - Conditional expectation
                    %     inverse.L      - Lead field (reduced)
                    %     inverse.R      - Re-referencing matrix
                    %     inverse.qC     - spatial  covariance
                    %     inverse.qV     - temporal correlations
                    %     inverse.T      - temporal subspace (discrete cosine set of x t temporal modes or frequencies)
                    %     inverse.U      - spatial  subspace
                    %     inverse.Is     - Indices of active dipoles
                    %     inverse.Nd     - number of dipoles
                    %     inverse.pst    - pers-stimulus time
                    %     inverse.dct    - frequency range
                    %     inverse.F      - log-evidence
                    %     inverse.R2     - variance accounted for (%)

                    P.trials=cond;
                    P.type=settings.Method;
                    P.smooth=settings.SourceSmoothness; % depends on mesh resolution. 0.6 (default) good for 4004 dipoles
                    P.Np=settings.SparsePriorsPerHemisphere; % See Friston et al submitted 07. Per hemisphere...so total=Np*3? Also depends on mesh resolution?
                    try
                        if isnumeric(settings.MaxSpatialModes); Nm=settings.MaxSpatialModes;
                        elseif ischar(settings.MaxSpatialModes); Nm=D.MaxFilter.HarmonicComponents;
                        else Nm=80; fprintf('\nUsing default maximum number of %g spatial modes\n',Nm)
                        end
                    catch Nm=80; fprintf('\nUsing default maximum number of %g spatial modes\n',Nm)
                    end
                    P.Nm=Nm; % (def:96) 86? See Nenonen 2007
                    P.Nr=[]; % allow more temporal modes?? def: 8
                    P.Na=Inf; % Want all data points for doing l-r subtraction
                    P.sdv = 4; % does this depend on sample rate? No: stdev of temporal autocorrelation in ms
                    P.lpf=-Inf; % remove lpf?? (1Hz by default) This is important. If>-Inf the mean across the window will be removed.
                    P.hpf=45; % (256Hz by default but max 48 for MEG??)
                    P.Han=0; % 0 to avoid source waveforms starting and ending at same value?
                    if tw==1 % include baseline period in 1st time window
                        P.woi=[-D.events.start*1000/D.Radc TimeWindows(tw,2)];
                    else P.woi=TimeWindows(tw,1:2);
                    end

                    %% create ID string for these inversion settings
                    temp=struct2cell(orderfields(rmfield(settings,{'Overwrite','ContrastType'})));
                    id='';
                    for c=2:length(temp);
                        id=[id regexprep(num2str(temp{c}),'\s*',',') '_'];
                    end
                    id=sprintf('%s%s_e%g_',id,mat2str(P.woi),condnum);
                    if dofuse; id=[id 'fusion']; end
                    id=regexprep(id,{'/','\'},{'',''});

                    status='unknown'; tempfm='';
                    %% check if this inversion needs to be done for this forward model, finding
                    %% correct index, or copying forward model to new one if necessary
                    for v=1:length(D.inv)
                        if strcmp(D.inv{v}.forward.method,forwardmodels{fm})
                            if ~isfield(D.inv{v},'comment') || isempty(D.inv{v}.comment)
                                status='todo'; break; % id field has not been
                                % created, so assume not yet inverted
                            end
                            D.inv{v}.comment=regexprep(D.inv{v}.comment,{'/','\'},{'',''});
                            try
                                D.inv{v}.inverse.id=regexprep(D.inv{v}.inverse.id,{'/','\'},{'',''});
                            catch
                            end
                            if isstruct(D.inv{v}.inverse) % save as sep file to stop D bloating
                                inverse=D.inv{v}.inverse;
                                fout=fullfile(D.path,'inversions',D.inv{v}.comment);
                                save(fout,'inverse','-MAT');
                                D.inv{v}.inverse=fullfile(D.path,'inversions',D.inv{v}.comment);
                            end
                            if strcmp(D.inv{v}.comment,id);
                                % inversion with these settings already done
                                if settings.Overwrite; status='todo';
                                else status='done';
                                end
                                break
                            elseif isempty(tempfm);
                                % this was inverted with different settings; store
                                % in case need to append new copy of forward model
                                % for this inversion, and keep looking
                                tempfm=D.inv{v}; tempfm.inverse='';
                                tempfm.comment=''; tempfn.contrast={};
                            end
                        end
                    end

                    if strcmp(status,'unknown')
                        v=v+1; D.inv{v}=tempfm;
                    end
                    condinvs(condnum)=v;

                    if ~strcmp(status,'done')
                        D.val=v;
                        %% do the inversion
                        %addpath /imaging/local/spm/spm5/cbu_updates
                        D.inv{D.val}.inverse=P;

                        colormap gray
                        fprintf('Forward model: %s, inverting with id: %s\n', ...
                            forwardmodels{fm},id)

                        % Inversion seems not to cope with bad channels. With trans
                        % to default in one step, MEG signal can get very large
                        % leading to channels being marked bad (perhaps
                        % unnecessarily?) Unmark them here!
                        D.channels.Bad=[];

                        if ~dofuse % invert just mags or grads alone, or both with prespecified weighting
                            if isempty(D.data) % load data if necessary
                                save(fullfile(D.path,D.fname),'D');
                                D=spm_eeg_ldata(fullfile(D.path,D.fname));
                            end
                            if settings.RotatingDipoles
                                D.inv{D.val}.forward.gainmat=D.inv{D.val}.forward.gainxyz;
                            end

                            %try
                            D=spm_eeg_invert(D); % my version would catch an unspecified field
                            %catch
                            %    aas_log(aap,1,'\nInversion failed. Check memory and that the leadfield matrix was constructed properly.\n')
                            %end
                            %D=spm_eeg_invert_dm_old(D); % preload file array to hopefully avoid _utExeption
                        else % fusion inversion
                            Dm.inv{1}.inverse=P;
                            Dg.inv{1}.inverse=P;
                            Dm.val=1; Dg.val=1;
                            if settings.RotatingDipoles
                                Dm.inv{1}.forward.gainmat=Dm.inv{1}.forward.gainxyz;
                                Dg.inv{1}.forward.gainmat=Dg.inv{1}.forward.gainxyz;
                            end
                            F = spm_eeg_invert_fuse({Dm,Dg}); %
                            D.inv{D.val}=F.inv{1};
                        end

                        %% The penultimate step is to specify a time-frequency window and compute the
                        % Conditional expectation of the RMS response.  A simple windowed average
                        % is a special case of this, where the frequency is zero.  In this context,
                        % the RMS is the same as the absolute value of the time-averaged repsone.
                        try D.inv{D.val}.contrast{1}; catch; D.inv{D.val}.contrast={}; end
                        % By default only one time/frequency contrast is estimated per inversion,
                        % so estimate each contrast on a copy then append to the inversion.
                        % Do this both for evoked and induced effects as requested.
                        S=D;
                        S.inv=S.inv(S.val); S.val=1;
                        for fw=1:size(FrequencyWindows,1)
                            for induced=0:1
                                if induced && strcmp(settings.ContrastType,'evoked'); continue; end
                                if ~induced && strcmp(settings.ContrastType,'induced'); continue; end
                                % any other value to run both

                                [pth nam]=fileparts(S.fname);
                                S.inv{1}.contrast='';
                                warning off all
                                S.inv{1}.contrast.woi =TimeWindows(tw,1:2); % peristimulus time (ms)
                                S.inv{1}.contrast.fboi=FrequencyWindows(fw,:); % frequency window (Hz)
                                if max(S.inv{1}.contrast.fboi)==0; S.inv{1}.contrast.fboi=[]; end;
                                warning on all

                                %%%% don't bother breaking up frequencies or doing induced power for short time
                                %%%% window; probably want to (impr|rem)ove this!
                                if range(S.inv{1}.contrast.woi)<100;
                                    if induced || ~isempty(S.inv{1}.contrast.fboi); continue; end
                                end
                                %%%%

                                if induced; type='induced'; else type='evoked'; end
                                % could also get evoked effects by running 'induced' on
                                % averaged data

                                cname=sprintf('%s_%sms_%sHz_%s', nam, ...
                                    regexprep(mat2str(S.inv{1}.contrast.woi),{' ','[',']'},{'-','',''}), ...
                                    regexprep(mat2str(S.inv{1}.contrast.fboi),{' ','[',']'},{'-','',''}), ...
                                    type(1:3));

                                % find contrast if it exists else add new one
                                done=0;
                                clear tempval
                                for dc=1:length(D.inv{D.val}.contrast)
                                    try D.inv{D.val}.contrast{dc};
                                        if ~isfield(D.inv{D.val}.contrast{dc},'cname'); break;
                                        elseif strcmp(cname,D.inv{D.val}.contrast{dc}.cname); 
                                            if evoked, tempval=dc; end % for combining with induced power later
                                            done=1; break;
                                        end
                                    catch; keyboard; break % should not happen?
                                    end
                                    if dc==length(D.inv{D.val}.contrast); dc=dc+1; end
                                end
                                if ~settings.Overwrite && done,
                                    fprintf('\nFound contrast window %g: %s, with %g named contrasts', ...
                                        dc,cname,length(D.inv{D.val}.contrast{dc}.names));
                                    continue;
                                end
                                S.inv{1}.contrast.type=type;
                                S.inv{1}.contrast.svd=settings.svd;
                                S = spm_eeg_inv_results_dm(S); % why does this change S.val?? (e.g. from 1 to 2)
                                % with spm_eeg_inv_results, if using a frequency band,
                                % the fourier set is truncated by svd.
                                % spm_eeg_inv_results_dm gives more options (not sure which
                                % is best!)...
                                % S.svd=0: turn it off
                                % S.svd=1: leave it as it is (to use eigenvectors)
                                % S.svd=2: do it, but then reweight eigenvectors by eigenvalues
                                % (also offers unrecommended rectangular time window,
                                % and loads data if necessary)

                                % copy to new field which I'll zscore across mesh
                                S.inv{1}.contrast.MW={single(zscore(S.inv{1}.contrast.GW{1}))};
                                % and do no scaling when writing volumes

                                % add contrast names
                                S.inv{1}.contrast.names=enames(condnum);
                                S.inv{1}.contrast.cname=cname;

                                if induced && exist('tempval','var'); % check if evoked has been done too and sum the power;
                                    % this measure of 'total power' is not
                                    % strictly correct (Friston '06) but I'm
                                    % curious...
                                    S.inv{1}.contrast.TW={single(zscore( ...
                                        S.inv{1}.contrast.GW{1} + D.inv{D.val}.contrast{tempval}.GW{1}))}
                                end
                                
                                % and append contrast to inversion
                                try D.inv{D.val}.contrast{dc}=S.inv{1}.contrast;
                                catch D.inv{D.val}.contrast{1}=S.inv{1}.contrast;
                                end
                                
                            end % now induced activity?
                        end % next frequency band
                        clear S

                    end % already done?

                    % save some stats
                    if isstruct(D.inv{v}.inverse); I.inverse=D.inv{v}.inverse;
                    else I=load(D.inv{v}.inverse,'-MAT');
                    end
                    tits={sprintf('R2:%s',id), sprintf('F:%s',id)};
                    vals={I.inverse.R2; I.inverse.F};
                    defs={'Percent variance explained by inversion.';...
                        'Log evidence of inversion (2nd step - MSP complexity ignored)'};
                    if min(P.woi)<0 % include mean SNR if there is baseline in inversion
                        stc=I.inverse.J{1}*I.inverse.T';
                        leng=fix(range(TimeWindows(tw,:)*D.Radc/1000));
                        start=fix(TimeWindows(tw,1)*D.Radc/1000)+D.events.start;
                        sig=mean(std(stc(:,start:start+leng),1,2)); % mean across vertices, of sample std across contrast window
                        noise=mean(std(stc(:,1:min(fix(0.1*D.Radc),leng)),1,2)); % std over min(100ms,same # samples from start of baseline)
                        tits=[tits, sprintf('SNR:%s',id)];
                        vals=[vals; sig/noise];
                        defs=[defs; 'Sample sd over contrast window/Sample sd over first 100ms of baseline, averaged over vertices'];
                    end
                    aas_emeg_savestats(D,tits,vals,defs);

                    if strcmp(status,'done');
                        fprintf('Forward model: %s, already inverted with id: %s\n', ...
                            forwardmodels{fm},id)
                        continue
                    end

                    %% write id to indicate completion with these settings, and save
                    D.data=[];
                    D.inv{v}.inverse.id=id;
                    D.inv{v}.comment=id;
                    fprintf('\nSaving...')
                    if isstruct(D.inv{v}.inverse) % save as sep file to stop D bloating
                        inverse=D.inv{v}.inverse;
                        fout=fullfile(D.path,'inversions',D.inv{v}.comment);
                        save(fout,'inverse','-MAT');
                        D.inv{v}.inverse=fullfile(D.path,'inversions',D.inv{v}.comment);
                    end
                    outfile=fullfile(D.path,D.fname);
                    save(outfile,'D');
                    copyfile(outfile,strrep(outfile,'.mat','.bck'),'f'); % backup

                    %% print summary of dipole activity for all events and inversions?
                    % if ~exist(fullfile(D.path,[nam '_sources.ps']),'file'),
                    % source_plots(D); end

                end % next event type
            end % next time window

            % add contrasts of unsigned localisations (coul)
            % and contrasts of signed localisations (cosl)
            % and average of signed localisations (aosl)?
            % Unsigned (zeroed) contrasts of signed localisations
            % (zucosl) seemed to give most interesting results, but
            % seem invalid because 2nd level test will not necessarily have zero
            % mean under H0.

            for c=1:length(Ccol)
                con=n(:,Ccol(c))'*n(:,Ecol);

                for fb=1:length(D.inv{v}.contrast) % v should be the last inverted condition for these settings
                    iid=regexprep(D.inv{condinvs(1)}.comment,'_e\d\d?_',sprintf('_c%g_',length(Ecol)+c));
                    cname=D.inv{condinvs(1)}.contrast{fb}.cname;
                    if strcmp('ind',cname(end-2:end)); schemes={'coip','cotp'}; % contrasts of induced and 'total' power
                    else schemes={'coep','coeJ'}; % contrasts of evoked power and evoked signed source estimate
                    end
                    for s=1:length(schemes)
                        status='unknown';
                        %% check if this contrast needs to be done for this
                        %% inversion type
                        for next_s=1:length(D.(schemes{s}))
                            if strcmp(D.(schemes{s}){next_s}.contrast{1}.cname,cname) ...
                                    && strcmp(D.(schemes{s}){next_s}.comment,iid);
                                % contrast with these settings already done
                                if settings.Overwrite; status='todo';
                                else status='done';
                                end
                                break
                            end
                        end
                        switch status
                            case 'unknown'
                                if isempty(next_s); next_s=1;
                                else next_s=next_s+1; % add new contrast
                                end
                            case 'done'
                                continue
                            case 'todo'
                                % overwrite current contrast
                        end

                        C=zeros(size(D.inv{v}.contrast{fb}.JW{1}));
                        for e=1:length(Ecol)
                            if strcmp(schemes{s},'coeJ')
                                C=C+con(e).*D.inv{condinvs(e)}.contrast{fb}.JW{1};
                            elseif strcmp(schemes{s},'coep')
                                C=C+con(e).*D.inv{condinvs(e)}.contrast{fb}.GW{1};
                            elseif strcmp(schemes{s},'cotp') && isfield(D.inv{condinvs(e)}.contrast{fb},'TW')
                                % note, this value of 'total power' not strictly appropriate (Friston 06) 
                                C=C+con(e).*D.inv{condinvs(e)}.contrast{fb}.TW{1};
                            end
                        end

                        D.(schemes{s}){next_s}.contrast{1}.MW{1}=zscore(C);
                        D.(schemes{s}){next_s}.contrast{1}.names{1}=enames{length(Ecol)+c};
                        D.(schemes{s}){next_s}.contrast{1}.cname=cname;
                        D.(schemes{s}){next_s}.comment=iid;
                    end % next contrast scheme
                end % next 1st level contrast
            end % next inversion
            fprintf('\nSaving...')
            save(fullfile(D.path,D.fname),'D');

        end % next forward model
    end % do fusion inversion too?
end % next file

