function [aap resp]=aamod_emeg_tf2vol(varargin)
% Average tf data over specified channels and write to freq x time vols
% Danny Mitchell 17/08/08

%% check task settings, subject, block etc
[aap subblock doit resp settings]=aa_emeg_checktasksettings(mfilename('fullpath'),varargin);
if ~doit; return; end

%% find files and decide whether to run task;
files=aas_emeg_findfiles(aap,settings.InputFilter,subblock);
if isempty(files); aas_log(aap,1,sprintf('\nFound no data! (Input filter is %s)\n',settings.InputFilter)); end

for f=1:length(files);

    %% load header
    try rehash; load(files{f});
    catch
        fprintf('\nRetrying in 60s in case of access conflict...')
        pause(60)
        try rehash; load(files{f});
        catch; aas_log(aap,1,'Failed! File may be corrupt? \n');
        end
    end
    fprintf('\nFile: %s...\n',D.fname)

    clear S
    S.Fname=files{f};
    S.overwrite=settings.Overwrite;

    %% average over specified channels
    % e.g. for channels from FDR significant region for main effect of
    % set-size for magnetometers from 0508 analysis:
    % settings.Channels={'MEG1631','MEG1841','MEG1911'};

    if iscell(settings.Channels)
        % get vector of channel indices from channel names
        % (all sensor types at that location)
        temp=[];
        for sn=1:length(settings.Channels)
            temp=[temp, find(strncmpi(D.channels.name,settings.Channels{sn}(1:6),6))];
        end
        changroups={temp};
    elseif ~isfield(settings,'Channels') || isempty(settings.Channels)
        %         % average all channels for each hemisphere
        [LHS, RHS]=getLRpairs2(D,'force_sideonly');
        %         % these are ordered anterior to posterior, but channels details are
        %         % stored alphabetically; don't think this matters here.
        %         [junk ind]=sort(D{k}.channels.name(LHS));
        %         LHS=LHS(ind);

        changroups={LHS,RHS};
        sufs={'-LHS','-RHS'}; % hyphen (not 'minus') confusing but consistent with suffix format for mags and grads
    else
        changroups={settings.Channels};
    end

    %% write volumes
    % accepts multiple files, but do one at a time as the channel indices will
    % be different if files contain different sensor types.
    for cg=1:length(changroups)
        fprintf('\nChannel set: %s:\n',sufs{cg})
        S.chans=changroups{cg};
        try S.suffix=sufs{cg}; catch; end
        fnames=spm_eeg_convertmat2ana3Dtf_dm(S);
    end

    %% smooth tf images
    try smooth=[str2num(settings.smooth) 0]; % 3rd dimension doesn't exist
    catch smooth=[];
    end
    if ~isempty(smooth)
        for tfi = 1:length(fnames)
            Pout          = fnames{tfi}; % don't bother with prefix
            if ~exist(Pout,'file') || settings.Overwrite
                fprintf('\nSmoothing %s',fnames{tfi})
                spm_smooth(spm_vol(fnames{tfi}),Pout,smooth);
            else fprintf('.')
            end
        end
    end

    if ~isfield(settings,'Channels') || isempty(settings.Channels)
        % calculate difference: left - right sensors, and average;
        for trialtype=1:length(fnames)
            Pin=char(strrep(fnames{trialtype},'RHS.img','LHS.img'),fnames{trialtype});
            func1='i1-i2';
            func2='(i1+i2)/2';
            Pout=strrep(fnames{trialtype},'RHS.img','L-R.img');
            if ~exist(Pout,'file') || settings.Overwrite
                fprintf('\nCalculating %s',Pout)
                spm_imcalc_ui(Pin,Pout,func1);
            end
            Pout=strrep(fnames{trialtype},'RHS.img','L+R.img');
            if ~exist(Pout,'file') || settings.Overwrite
                fprintf('\nCalculating %s',Pout)
                spm_imcalc_ui(Pin,Pout,func2);
            end
        end
    end

end % next file

return
