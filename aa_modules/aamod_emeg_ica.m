function [aap resp]=aamod_emeg_ica(varargin)
% Calculate independent components
% Adapted from Jason's _example_pathway.m
%
% Danny Mitchell 04/04/08

%% check task settings, subject, block etc
[aap subblock doit resp settings]=aa_emeg_checktasksettings(mfilename('fullpath'),varargin);
if ~doit; return; end

%% find files and decide whether to run task;
files=aas_emeg_findfiles(aap,settings.InputFilter,subblock);
if isempty(files); aas_log(aap,1,sprintf('\nFound no data! (Input filter is %s)\n',settings.InputFilter)); end

%% add paths for eeglab
addpath /imaging/local/spm_eeglab/
addpath /imaging/local/eeglab

%%
for f=1:length(files);

%% Run ICA (on continuous or epoched data):
    % This could take some time, 20 minutes - ? hours for continuous data?!
    % If you hit memory problems, adjust S.samppct down.

    [pth nam ext]=fileparts(files{f});
    logfile=fullfile(pth,['i' nam '.log']);
    i_file=fullfile(pth,['i' nam ext]);
    if ~exist(strrep(logfile,'.log','.done'),'file') || settings.Overwrite
        D=spm_eeg_ldata(files{f});
        clear S
        S.D=files{f};
        S.samppct=settings.DataFraction;
        S.newpfx='i';
         %try ncomps=D.MaxFilter.HarmonicComponents; 
         %catch % restrict components to minimise memory problems
             %
             % EEGLAB tutorial recommends > 30 samples per weight point
             % let's use pca to ensure 50 samples per weight point;
             % setting 'pca' to 'auto' will do this.
             % (For 15min continuous data, dowsampled x4, will be ~64;
             % 64 is number of inner components used by maxfilter with
             % lin=8, lout=3)
         %end
         if ischar(settings.pca); pca='auto'; else pca=settings.pca; end
        S.args={'extended',1,'pca',pca,'maxsteps',1000};
        S.excrej=1; % assuming we want to exlude rejected trials
        S.excbad=1; % assuming we want to exlude rejected channels
        if exist(logfile,'file'); delete(logfile); end
        diary(logfile)
        D=spm_eeglab_runica_dm(S); 
        % Don't duplicate data file; select input data differently; use
        % binica rather than runica.
        diary off
        movefile(logfile,strrep(logfile,'.log','.done'));
    else load(i_file);
    end
    clear S acts

%% Get ICA activations:
    Ai_file=fullfile(pth,['Ai' nam ext]);
    if ~exist(Ai_file,'file') || settings.Overwrite
        S.D=i_file;
        S.samppct=1;
        % if adjust samppct, also need to adjust eogdata accordingly:
        % e.g., eogdata=eogdata(:,1:ceil(D.Nsamples*S.samppct));
        S.newfname=['Ai' nam ext];
        S.excrej=1; % assuming we want to exlude rejected trials
        fprintf('\nCalculating ICs...');
        spm_eeglab_icaact_dm(S); % trying to save memory and deal with file array problems if rejecting epochs?
    end
    
    % %%% Plot activations and EOG (child window):
    % % Build child-window structure:
    % eogdata=D.data([D.channels.veog D.channels.heog],:);
    % C.D=D.fname;
    % C.data=eogdata;
    % C.args={'winlength',10,'spacing',500};
    % chfig=spm_eeglab_eegplot(C);% Plot child window:
    % % Build main-window structure:
    % S.D=D.fname;
    % S.data=acts;  % plots ICactivations instead of data!
    % S.args={'winlength',10,'spacing',10,'dispchans',12,'children',chfig};
    % mainfig=spm_eeglab_eegplot(S);% Plot main window:

end % next file

return
