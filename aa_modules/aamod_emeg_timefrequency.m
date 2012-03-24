function [aap resp]=aamod_emeg_timefrequency(varargin)
% Compute time-frequency decomposition of power or magnitdude, and
% phase-locking-value at sensor level. (Quite slow.)
% Can run on averaged file to get time-frequency plot of evoked power
%
% Modified spm_eeg_tf to
% - show progress bar,
% - convolve all channels at once for speed
% - add optional extra label to output file name (although no longer use
% this; was to mark mags vs grds or induced vs evoked)
% - average over events as well to avoid memory problems
%
% For each event, first convolves all frequencies and channels...
% Then subtracts mean baseline power per epoch if requested...
% Then collapses over channels if requested.
%
% How sensible is it to run this on filtered data?
% Could do this at source level??
% What does PLV mean for averaged data?
%
% To avoid edge effects, baseline and analysed timepoints should be at
% least half the wavelet width away from the edge of the epoch
% I.e. > N/2 * 1000/Freq., where N is Morlet order.
% E.g. for 4Hz at N=6: > 6/2 * 1000/4
%                      > 750ms
% Also, any baseline period should be at least a full cycle wide
% e.g. for 4Hz, > 125ms
% Therefore, effectively, in this example a 875ms pre-trigger period would 
% be desireable (excluding further considerations that baseline window 
% might be best kept similarly seperate from signal of interest). 
% See e.g. Roach 2008.
%
% Danny Mitchell 03/03/08

%% check task settings, subject, block etc
[aap subblock doit resp settings]=aa_emeg_checktasksettings(mfilename('fullpath'),varargin);
if ~doit; return; end

%% find files and decide whether to run task;
files=aas_emeg_findfiles(aap,settings.InputFilter,subblock);
if isempty(files); aas_log(aap,1,sprintf('\nFound no data! (Input filter is %s)\n',settings.InputFilter)); end

%% run task for each file

S=settings;
S.label='';
S.frequencies=eval(S.frequencies); % this is entered as string e.g. '4:40'
S.pow=(S.pow==1); % this is entered as 2 for magnitude instead of 0
% frequencies   - vector; Hz; non-linear steps don't seem to be plotted correctly;
%                   lower bound must be >0
% rm_baseline	- subtract baseline power (1/0) yes/no
% Mfactor       - Morlet wavelet factor (can not be accessed by GUI)
%                   default is 7; Rik used 5; Wikipedia suggests>5;
%                   trades frequency vs time resolution with lower orders favouring time
% pow           - 1 = power, 0 = magnitude
% collchans     - 1 = collapse across channels, 0 = do not

for f=1:size(files,1);

    % get data file and channels and baseline period as necessary
    S.D=files{f};
    [pth,nam]=fileparts(S.D);
    localdoneflag=fullfile(pth, ['mt2' S.label '_' nam '.done']);
    
    if ~exist(localdoneflag,'file') || settings.Overwrite==1;
        if exist(localdoneflag,'file'); delete(localdoneflag); end;
        load(S.D);
        if isinf(settings.channels); % select all emeg channels
            S.channels=D.channels.eeg; 
        end
        if S.rm_baseline
            if ischar(settings.Sbaseline); 
                S.Sbaseline(1)=1;
                S.Sbaseline(2)=D.events.start;
            end
        end
    
        fprintf('\nProcessing file: %s',files{f})
        clear D pth nam ext
        D=spm_eeg_tf_dm3(S);
        % Modified spm_eeg_tf to
        % - show progress bar,
        % - convolve all channels at once for speed
        % - add extra label to output file name (currently unused)
        % - v3: do the averaging in here to avoid memory problems
        
        diary(localdoneflag); diary off;
    end

    %%% THIS SHOULD BE MOVED ABOVE DIARY(LOCALDONEFLAG)
        % add mask to exclude regions where there might be edge effects
        temp=fullfile(pth,['mt1_' nam]);
        load(temp);
        
        buffer=round(S.Mfactor/2 .* 1000./S.frequencies * D.Radc/1000); % in samples
        
         D.tf.mask=zeros(D.Nfrequencies,D.Nsamples);
        for cf=1:D.Nfrequencies
            D.tf.mask(cf,buffer(cf):end-buffer(cf))=1;
        end
        save(fullfile(D.path,D.fname),'D');
        
        load(strrep(temp,'mt1','mt2'))
        
         D.tf.mask=zeros(D.Nfrequencies,D.Nsamples);
        for cf=1:D.Nfrequencies
            D.tf.mask(cf,buffer(cf):end-buffer(cf))=1;
        end
        save(fullfile(D.path,D.fname),'D');
    %%%
    
    fprintf('.')
end % next file

return