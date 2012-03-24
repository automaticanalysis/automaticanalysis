function [aap resp]=aamod_emeg_contrasts(varargin)
% Run spm_eeg_average first to apply any artefact rejection/weighting and
% get ERPs/ERFs.
% Then run spm_eeg_weight_epochs to apply any contrasts.
% 
% Danny Mitchell 03/03/08

%% check task settings, subject, block etc
[aap subblock doit resp settings]=aa_emeg_checktasksettings(mfilename('fullpath'),varargin);

%% CHECK REQUIREMENTS: load spreadsheet containing epoching instructions
epochfile=spm_select('FPList',fullfile(aas_getsesspath(aap,subblock{1},subblock{2}),'events'),'(E|e)pochs.*.xls');
if exist(epochfile,'file')~=2; epochfile=spm_select('FPList',fullfile(aas_getsubjpath(aap,subblock{1}),'events'),'(E|e)pochs.*\.xls'); end
    if exist(epochfile,'file')~=2; epochfile=spm_select('FPList',fullfile(aap.acq_details.root,'events'),'(E|e)pochs.*.xls'); end
    if exist(epochfile,'file')~=2; error('aa:EpochFileError', '\n Failed to find Epochs.*.xls, required for epoching continuous data. \n This should be placed in "events" folder at either study, subject, or session levels, with lower levels overriding higher ones. \n See examples/Epochs.xls for an example of the required format. \n'); end
    warning off all; [Numeric,Txt]=xlsread(epochfile); warning on all
    
if ~doit; return; end

%% find files and decide whether to run task;
files=aas_emeg_findfiles(aap,settings.InputFilter,subblock);
if isempty(files); aas_log(aap,1,sprintf(' Found no data! Using input filter "%s"',settings.InputFilter)); return; end

%% do it

for f=1:length(files);

    load(files{f});

    [pth name ext]=fileparts(files{f});
    c=Numeric(ismember(Numeric(:,1),D.events.types),:);
    fc=find(any(c(:,3:end)))+2;

    if ~isfield(D.events,'names') % add event names if not there already
        fe=sum(c==0)==size(c,1)-1 & sum(c)==1;
        D.events.names=Txt(1,fe);
        save(fullfile(D.path,D.fname),'D');
    end
    S.D=files{f};
    outputfile=fullfile(pth,['m' name ext]);

    if ~exist(outputfile,'file') || settings.Overwrite==1 % || 1==1;
        %addpath /imaging/dm01/MEG/aaMEG % reconfirm my version of spm_eeg_ldata is found, for spm_eeg_average

        fprintf('\nAveraging %s...',outputfile);
        % Average events including artefact weighting/rejection.
        % This can be quite slow.
        D=spm_eeg_average_dm(S); % Written to file with m prefix.
        S.c=c(:,fc)'; % create contrast matrix for averaged data

        % Always average first so artefact rejection/correction is
        % applied!
        %             % create contrast matrix for epoched data
        %             S.c=zeros(length(fc),D.Nevents);
        %             for con=1:length(fc)
        %                 for ev=1:length(D.events.types)
        %                     S.c(con,D.events.code==D.events.types(ev))=Numeric(Numeric(:,1)==D.events.types(ev),fc(con));
        %                 end
        %             end
        %             % add vector of ones to repl field to indicate that each epoch is a single trial
        %             D.events.repl=ones(1,D.Nevents);
        %             save(fullfile(D.path, D.fname), 'D');

        % apply the contrast
        fprintf('\nApplying contrasts...');
        S.D=fullfile(D.path,D.fname);
        S.events.WeightAve=0; % hmm
        D = spm_eeg_weight_epochs_dm(S); % this could average epochs as well if not done already
        % but doesn't apply robust weights or artefact rejections. Adds
        % another m prefix. Remove this now...
        cmd=['mv ' fullfile(D.path,D.fname) ' ' fullfile(D.path,D.fname(2:end))]; unix(cmd);
        cmd=['mv ' fullfile(D.path,D.fnamedat) ' ' fullfile(D.path,D.fnamedat(2:end))]; unix(cmd);
        D.fname=D.fname(2:end);
        D.fnamedat=D.fnamedat(2:end);
  
        % add event & contrast names from spreadsheet
        D.events.names=Txt(1,fc);
        save(fullfile(D.path,D.fname),'D');   
        
        % save SNR to mat file
        if ~isfield(D,'data') || isempty(D.data);
            D=spm_eeg_ldata(fullfile(D.path,D.fname));
        end
        for e=1:size(D.data,3)
            sig=mean(std(D.data(D.channels.eeg,D.events.start:end,e),1,2)); % mean across channels, of sample std over dim 2
            noise=mean(std(D.data(D.channels.eeg,1:0.1*D.Radc,e),1,2)); % use 1st 100ms of baseline
            snr=sig/noise;
            tits={sprintf('Signal%g',e); sprintf('Noise%g',e); sprintf('SNR%g',e)};
            vals={sig;noise;snr};
            defs={'Sample standard deviation over first 100ms of baseline, averaged over sensors and stimulus types.';...
                'Sample standard deviation over all post-trigger timepoints, averaged over sensors and stimulus types.';...
                'Empirical SNR estimate, similar to Henson paper on EMEG fusion.'};
            if strcmp(D.modality,'MEG')
                Msig=mean(std(D.data(1:102,D.events.start:end,e),1,2)); % mean across channels, of sample std over dim 2
                Mnoise=mean(std(D.data(1:102,1:0.1*D.Radc,e),1,2)); % use 1st 100ms of baseline
                Msnr=Msig/Mnoise;
                Gsig=mean(std(D.data(103:end,D.events.start:end,e),1,2)); % mean across channels, of sample std over dim 2
                Gnoise=mean(std(D.data(103:end,1:0.1*D.Radc,e),1,2)); % use 1st 100ms of baseline
                Gsnr=Gsig/Gnoise;
                tits=[tits; {sprintf('SNRm%g',e); sprintf('SNRg%g',e)}];
                vals=[vals;{Msnr; Gsnr}];
                defs=[defs; {'SNR for magnetometers only.'; 'SNR for gradiometers only.'}];
            end

            aas_emeg_savestats(D,tits,vals,defs,settings.Overwrite);
        end
    end

    % plot timecourses and plot topographies moved to their own modules

    fprintf('.')
end

return