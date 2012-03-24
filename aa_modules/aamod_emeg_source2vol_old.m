function [aap resp]=aamod_emeg_source2vol(varargin)
% Uses spm_eeg_inv_results to compute conditional expectation of the RMS 
% response across time/frequency windows specified in events spreadsheet,
% then writes results to volume with 2D and 3D smoothing.
%
% Danny Mitchell 02/04/08

%% check task settings, subject, block etc
[aap subblock doit resp settings]=aa_emeg_checktasksettings(mfilename('fullpath'),varargin);

%% CHECK REQUIREMENTS: load spreadsheet containing epoching instructions
epochfile=spm_select('FPList',fullfile(aas_getsesspath(aap,varargin{3},varargin{4}),'events'),'(E|e)pochs.*.xls');
if exist(epochfile,'file')~=2; epochfile=spm_select('FPList',fullfile(aas_getsubjpath(aap,varargin{3}),'events'),'(E|e)pochs.*\.xls'); end
if exist(epochfile,'file')~=2; epochfile=spm_select('FPList',fullfile(aap.acq_details.root,'events'),'(E|e)pochs.*.xls'); end
if exist(epochfile,'file')~=2; error('aa:EpochFileError', '\n Failed to find Epochs.*.xls, required for epoching continuous data. \n This should be placed in "events" folder at either study, subject, or session levels, with lower levels overriding higher ones. \n See examples/Epochs.xls for an example of the required format. \n'); end

if ~doit; return; end

%% find files and decide whether to run task;
files=aas_emeg_findfiles(aap,settings.InputFilter,subblock);
if isempty(files); aas_log(aap,1,sprintf('\nFound no data! (Input filter is %s)\n',settings.InputFilter)); end

%% load event specifications from spreadsheet
warning off all; [Numeric,Txt]=xlsread(epochfile); warning on all

for f=1:length(files);

%% skip EEG or CDA files
    if ~isempty(regexp(files{f},'-eeg','ONCE')); continue; end
    if ~isempty(regexp(files{f},'_CDA_','ONCE')); continue; end
    fprintf('.')
    fprintf('\nFile %s:',files{f});
%% load MEG header
    try rehash; load(files{f});
    catch
        fprintf('\nFailed to load file. Will try again in 60s in case there is an access conflict...\n');
        pause(60)
        try rehash; load(files{f},'-MAT');
        catch; 
            %aas_log(aap,1,'\n File may be corrupt? \n');
            fprintf('\nDebug %s\n',mfilename); keyboard
        end
    end

%% load relevant events/contrasts
    try
        % worked when events were labeled from contrast header row
        % temptext=regexprep(Txt(1,:),sprintf('%s|',D.events.names{:}),'match');
        % col=strcmp(temptext,'match');
        
        % works now that events might be labeled from event name column
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
    Ecol=find(size(n,1)-sum(n==0,1)==1); % relevant columns containing events
    Ccol=find(round(sum(n))==0); % relevant columns containing zeroed contrasts
    
    % for each inversion
    for v=1:length(D.inv);
        D.val=v;
        fprintf('\n - Inversion %g:',v);
%% The penultimate step is to specify a time-frequency window and compute the
        % Conditional expectation of the RMS response.  A simple windowed average
        % is a special case of this, where the frequency is zero.  In this context,
        % the RMS is the same as the absolute value of the time-averaged repsone.
        try D.inv{D.val}.contrast{1}; catch; D.inv{D.val}.contrast={}; end
        % By default only one time/frequency contrast is estimated per inversion,
        % so estimate each contrast on a copy then append to the inversion.
        % Do this both for evoked and induced effects.
        S=D;
        for tw=1:size(TimeWindows,1)
            for fw=1:size(FrequencyWindows,1)
                for induced=0:1
                    [pth nam]=fileparts(S.fname);
                    S.inv{S.val}.contrast='';
                    warning off all
                    S.inv{S.val}.contrast.woi =TimeWindows(tw,1:2); % peristimulus time (ms)
                    S.inv{S.val}.contrast.fboi=FrequencyWindows(fw,:); % frequency window (Hz)
                    if max(S.inv{S.val}.contrast.fboi)==0; S.inv{S.val}.contrast.fboi=[]; end;
                    warning on all

                    %%%% don't bother breaking up frequencies for short time
                    %%%% window; probably want to (impr|rem)ove this!
                    if ~isempty(S.inv{S.val}.contrast.fboi) && range(S.inv{S.val}.contrast.woi)<100; continue; end
                    %%%%

                    if induced; type='induced'; else type='evoked'; end
                    % could also get evoked effects by running 'induced' on
                    % averaged data

                    cname=sprintf('%s_%sms_%sHz_%s', nam, ...
                        regexprep(mat2str(S.inv{S.val}.contrast.woi),{' ','[',']'},{'-','',''}), ...
                        regexprep(mat2str(S.inv{S.val}.contrast.fboi),{' ','[',']'},{'-','',''}), ...
                        type(1:3));

                    % find contrast if it exists else add new one
                    done=0;
                    for dc=1:length(D.inv{D.val}.contrast)
                        try D.inv{D.val}.contrast{dc};
                            if ~isfield(D.inv{D.val}.contrast{dc},'cname'); break;
                            elseif strcmp(cname,D.inv{D.val}.contrast{dc}.cname); done=1; break;
                            end
                        catch; keyboard; break % should not happen?
                        end
                        if dc==length(D.inv{D.val}.contrast); dc=dc+1; end
                    end
                    if ~settings.Overwrite && done,
                        fprintf('\nFound contrast window %g: %s',dc,cname);continue;
                    end
                    S.inv{S.val}.contrast.type=type;
                    S.inv{S.val}.contrast.svd=settings.svd;
                    S = spm_eeg_inv_results_dm(S);
                    % with spm_eeg_inv_results, if using a frequency band,
                    % the fourier set is truncated by svd.
                    % spm_eeg_inv_results_dm gives more options (not sure which
                    % is best!)...
                    % S.svd=0: turn it off
                    % S.svd=1: leave it as it is (to use eigenvectors)
                    % S.svd=2: do it, but then reweight eigenvectors by eigenvalues
                    % (also offers unrecommended rectangular time window, 
                    % and loads data if necessary)

                    S.inv{S.val}.contrast.MW=S.inv{S.val}.contrast.JW;
                    % add power of any zeroed contrasts of signed localisations
                    for c=1:length(Ccol)
                        con=n(:,Ccol(c))'*n(:,Ecol);
                        zCoSL=zeros(size(S.inv{S.val}.contrast.JW{1}));
                        for e=1:length(Ecol)
                            zCoSL=zCoSL+con(e).*S.inv{S.val}.contrast.JW{e};
                        end
                        S.inv{S.val}.contrast.MW{end+1}=zCoSL.^2;
                    end

                    % zscore power to zero mean and unity standard deviation
                    % across mesh
                    for c=1:length(S.inv{S.val}.contrast.MW)                  
                        S.inv{S.val}.contrast.MW{c}=zscore(S.inv{S.val}.contrast.MW{c});
                    end

                    % add contrast names
                    S.inv{S.val}.contrast.names=Txt(1,col);
                    S.inv{S.val}.contrast.names=[ ...
                        regexprep(S.inv{S.val}.contrast.names(Ecol),'.*','$0(ZuL)') ...
                        regexprep(S.inv{S.val}.contrast.names(Ccol),'.*','$0(ZuCosL)')];
                    S.inv{S.val}.contrast.cname=cname;

                    % and append contrast to inversion
                    try D.inv{D.val}.contrast{dc}=S.inv{S.val}.contrast;
                    catch D.inv{D.val}.contrast{1}=S.inv{S.val}.contrast;
                    end

                end % evoked then induced
            end % next frequency window
        end % next time window
                    
        save(fullfile(D.path,D.fname),'D'); % is this best time to save? Can
        % be slow, so don't want to do it too often.
        
        %% Finally, write the smoothed contrasts to an image in normalised voxel space.
        for ctw=1:length(D.inv{D.val}.contrast)
            for induced=0:1
                if induced; type='induced'; else type='evoked'; end
                try cname=D.inv{D.val}.contrast{ctw}.cname;
                catch
                    [pth nam]=fileparts(S.fname);
                    cname=sprintf('%s_%sms_%sHz_%s', nam, ...
                        regexprep(mat2str(D.inv{D.val}.contrast{ctw}.woi),{' ','[',']'},{'-','',''}), ...
                        regexprep(mat2str(D.inv{D.val}.contrast{ctw}.fboi),{' ','[',']'},{'-','',''}), ...
                        type(1:3));
                end
                lastfile=fullfile(D.path,sprintf('sw_%s_%s_%02g.nii',cname,D.inv{D.val}.inverse.type,length(D.inv{D.val}.contrast{ctw}.MW)));
                if ~exist(lastfile,'file')|| settings.Overwrite==1
                    S=D;
                    S.inv{D.val}.contrast=D.inv{D.val}.contrast{ctw};
                    S.inv{D.val}.contrast.SmoothingSteps2D=settings.SmoothingSteps2D; 
                    S.inv{D.val}.contrast.SmoothingMethod=settings.Smoothing3D; 
                    % 'fixed' to apply smoothing as normal, or 'target' to 
                    % estimate mean smoothness and add extra smoothing to aim 
                    % for target final smoothness
                    S.inv{D.val}.contrast.smooth  = settings.Smoothing3DFWHM; % gaussian FWHM (mm; default 8? (not 12?); Rik used 16; follows graph lapacian smoothing on mesh)
                    S.inv{D.val}.contrast.display = 0;
                    S.inv{D.val}.contrast.overwrite = settings.Overwrite;
                    S = spm_eeg_inv_Mesh2Voxels_dm3(S);
                    % above uses MW but with no extra scaling;
                    % add contrast name to descrip;
                    % add t/f window and inverion type to file name
                    % option not to overwrite
                    % try to allow adjustable 3D smoothing for target smoothness
                    % option to set number of 2d smoothing steps
                    D.inv{D.val}.contrast{ctw}=S.inv{D.val}.contrast;

                    D.data=[]; save(fullfile(D.path,D.fname),'D');
                end
            end % do evoked then induced
        end % next time window
    end % next inversion

    save(fullfile(D.path,D.fname),'D'); % is this best time to save? Can
    % be slow, so don't want to do it too often.
        
end % next file

return
