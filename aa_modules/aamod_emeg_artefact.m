function [aap resp]=aamod_emeg_artefact(varargin)
% Apply robust averaging and/or threshold-based trial rejection.
%
% (The wtrials field from the spm manual doesn't seem to exist ???)
% Save weights to seperate file to stop header file bloating (then need to load with
% local copy of spm_eeg_ldata to memory map the file)
%
% Recommendations from James Kilner 26/11/07:
% Similar blocks should definitely be concatenated prior to robust averaging
% Always use the default weighting function offset
% Data should really be filtered with same filter settings as originally
% used rather than gaussian smoothing; might be better to use a small
% window then filter afterwards as desired. But floor(smoothing/2000*Radc)
% must be greater than 0.
%
% Danny Mitchell 03/03/08

%% check task settings, subject, block etc
[aap subblock doit resp settings]=aa_emeg_checktasksettings(mfilename('fullpath'),varargin);
if ~doit; return; end

%% find files and decide whether to run task;
files=aas_emeg_findfiles(aap,settings.InputFilter,subblock);
if isempty(files); aas_log(aap,1,' Found no data! '); return; end

%% run task for each file

for f=1:length(files);
    clear D S
    S.artefact=settings;
    S.D=files{f};
    S.thresholds.External_list=1;
    S.thresholds.out_list=[]; % trials to definitely reject
    S.thresholds.in_list=[]; % trials to definitely keep
        
%% prepare output file
    [pth name ext]=fileparts(files{f});
    if settings.RobustAveraging
        outputfile=fullfile(pth,['r' name ext]);
    else outputfile=fullfile(pth,['a' name ext]);
    end

    if ~exist(outputfile,'file') || settings.Overwrite==1;

%% plot trials rejected by EOG threshold, with optional automatic choice of rejection threshold
        D=spm_eeg_ldata(S.D);
        D.data=D.data(:,:,:);
        warning off all;
        Fh=spm_figure('GetWin','Graphics'); spm_figure('Clear',Fh);
        warning on all
        [S.artefact.thresholdVEOG veogbadlist]=showthresh(D,'VEOG',settings.thresholdVEOG,settings.thresholdVEOGmean);
        [S.artefact.thresholdHEOG heogbadlist]=showthresh(D,'HEOG',settings.thresholdHEOG,settings.thresholdHEOGmean);
        outplot=fullfile(pth,'figures',['a' name '_Rejection.png']);
        drawnow
        if D.channels.veog<D.channels.heog
            EOGthreshes=[S.artefact.thresholdVEOG S.artefact.thresholdHEOG];
        else
            EOGthreshes=[S.artefact.thresholdHEOG S.artefact.thresholdVEOG];
        end

%% prepare remaining parameters
           
        S.thresholds.out_list=union(S.thresholds.out_list,veogbadlist); 
        S.thresholds.out_list=union(S.thresholds.out_list,heogbadlist); % trials to definitely reject     
        
        % settings for robust averaging
        if settings.RobustAveraging==1
            S.artefact.weighted=1;
            S.artefact.weightingfunction=settings.WeightingFunction;
            S.artefact.smoothing=settings.Smoothing; % gaussian FWHM (over samples?)
        else S.artefact.weighted=0;
        end

        % thresholds
        if all(isinf([S.artefact.thresholdVEOG S.artefact.thresholdHEOG S.artefact.thresholdEEG S.artefact.thresholdMag S.artefact.thresholdGrad]))
            S.thresholds.Check_Threshold=0;
        else
            S.thresholds.Check_Threshold=1; % whether or not to reject by threshold
            if isempty(regexp(S.D,'-eeg','ONCE'))
                % could treat mags and grad the same (I think they should
                % have been weighted by now) but noise may well be
                % different.
                % max typical signal seems to be about 3000-5000 from a cursory glance
                S.thresholds.threshold=[ones(1,102)*S.artefact.thresholdMag ...
                    ones(1,204)*S.artefact.thresholdGrad ...
                    EOGthreshes];
            else
                if ~exist('D','var'); load(S.D); end
                S.thresholds.threshold=[ones(1,length(D.channels.eeg))*S.artefact.thresholdEEG ...
                    EOGthreshes];
            end
        end

%% do it and print diagnostic including rejected channels and trials
        %D=spm_eeg_artefact_dm(S); % just some small tweaks to hopefully speed up RA slightly
        S.usescaling=0; % for original import scheme; if using new scheme might need to change this
        D=spm_eeg_artefact(S);
        
        allrejected=find(D.events.reject);       
        rejected={};
        temp=allrejected;
        while length(temp)>10
            rejected=[rejected; num2str(temp(1:10))];
            temp(1:10)=[];
        end
        rejected=[rejected; num2str(temp)];
        
        % print #, proportion and list of rejected channels and events
        anntop=annotation(Fh,'textbox',[0 0.88 1 .12],'String', ...
            [{outplot;sprintf('Date:...%s',datestr(now,0)); ...
            sprintf('%g bad channels: %s',length(D.channels.Bad),num2str(D.channels.Bad)); ...
            sprintf('%g bad trials (%.0f%%, leaving %g):',length(allrejected), length(allrejected)/D.Nevents*100, D.Nevents-length(allrejected))};...
            rejected], ...
            'edge','none','fontsize',8);

        % find # of good trials per event
        for e=1:D.events.Ntypes
            goodperevent(e)=sum(~D.events.reject(D.events.code==D.events.types(e)));
        end

        % equalize trial number per event by random subsampling?
        norm=repmat({''},[D.events.Ntypes,1]);
        if settings.EquateEvents
            [least ind]=min(goodperevent);
            for e=1:D.events.Ntypes
                if e~=ind
                    numtoremove=goodperevent(e)-least;
                    good=find(D.events.code==D.events.types(e) & ~D.events.reject);
                    rand('state',sum(100*clock));
                    randgood=good(randperm(goodperevent(e))); %randomize
                    D.events.reject(randgood(1:numtoremove))=1;
                end
            end
            norm{ind}='_(all reduced to this)';
        end

        try eventnames=D.events.names(:);
        catch eventnames=cellstr(num2str(D.events.types(:)));
        end
        annbot=annotation(Fh,'textbox',[0 0 1 .5],'String', ...
            ['Good trials per event:',;strcat(eventnames,':_',cellstr(num2str(goodperevent(:))),norm)], ...
            'edge','none','fontsize',8);

        % print diagnostic
        drawnow
        print(Fh,'-dpng','-r82','-zbuffer',outplot);
        delete([anntop annbot]);

    else
        load(outputfile);
    end
    
%% memory map any weights from robust averaging
    if settings.RobustAveraging==1 && isfield(D,'weights')
        if ~isempty(D.weights)
            if isobject(D.weights); D.weights=[];
            else
                % write weights data to seperate file. This will then be memory
                % mapped by local copy of spm_eeg_ldata
                fprintf('Writing estimated weights to file...')
                D.weightsfname=strrep(D.fnamedat,'.dat','_weights.dat');
                D.weights=D.weights(:,:);
                fpd = fopen(fullfile(D.path, D.weightsfname), 'w');
                fwrite(fpd, D.weights, 'float32');
                fclose(fpd);
                D.weights=[];
            end
            save(fullfile(D.path,D.fname),'D');
        end

%% print diagnostic for robust averaging
        [pth nam]=fileparts(fullfile(D.path,D.fname));
        outplot=fullfile(pth,[nam '_RobAve.png']);
        if ~exist(outplot,'file') || settings.Overwrite==1;
            % prepare figure
            pos=0;
            warning off all;
            Fh=spm_figure('GetWin','Graphics'); spm_figure('Clear',Fh);
            warning on all
            % preload data from file array to avoid crashes
            addpath /imaging/dm01/MEG/aaMEG/
            D=spm_eeg_ldata(fullfile(D.path,D.fname));
            data=single(D.data(:,:,:));
            weights=reshape(D.weights(:,:),D.Nchannels,D.Nsamples,D.Nevents);
            for e=D.events.types
                pos=pos+1;
                % get one event type at time
                edata=data(:,:,D.events.code==e);
                eweights=weights(:,:,D.events.code==e);
                ewd=edata.*eweights;
                Nreps=sum(D.events.code==e);
                % find channel with max power over all samples and trials
                chanpow=sum(sum(edata.^2,3),2);
                [sorted ind]=sort(chanpow,'descend');
                cdata=spm_vec(squeeze(edata(ind(1),:,:)));
                cweights=spm_vec(squeeze(eweights(ind(1),:,:)));
                t=repmat(uint16(1:D.Nsamples)',[Nreps 1]);
                % plot each time point for each trial, colored by weight
                % subplot(ceil(D.events.Ntypes.^0.5),floor(D.events.Ntypes.^0.5),pos)
                subplot(D.events.Ntypes,1,pos)
                scatter(t,cdata,[],cweights,'.')
                colorbar
                axis tight
                % add mean across trials and robust average
                hold on
                plot(squeeze(mean(edata(ind(1),:,:),3)),'k','linewidth',2)
                plot(squeeze(sum(ewd(ind(1),:,:),3)./sum(eweights(ind(1),:,:),3)),'y')
                try ylabel(['Event:' D.events.names{pos}])
                catch ylabel(['Event trigger:' num2str(e)])
                end
                if pos==1; title('Example channels have max power over all samples and trials'); end;
                drawnow
                clear cdata cweights edata eweights ewd chanpow t
            end
            xlabel('Samples')
            legend({sprintf('Data from %g trials',Nreps),'Mean','Robust average'})
            colormap jet
            % print as bitmap
            fig=getframe(Fh);
            subplot('position',[0 0 1 1])
            image(fig.cdata); axis off image
            print('-dpng','-r82','-zbuffer',outplot)
        end
    end

    %% Rereference EEG without bad channels
    if strcmp(D.modality,'EEG') % && ~isempty(D.channels.Bad); ...might want to reference anyway
        fprintf('\nRereferencing EEG to mean of good channels...')
        if ~isfield(D,'data') || isempty(D.data)
            D=spm_eeg_ldata(outputfile);
        end
        meg_reavg_ref(D);
    end

end % next file

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [thresh badlist]=showthresh(D,whichchan,thresh,threshmean)
% e is data from whichchan
% whichchan is 'VEOG' or 'HEOG' or 'EEG'
% thresh is fixed value or 'auto'
if strcmp(whichchan,'VEOG'); 
    start=0; 
    e=squeeze(D.data(D.channels.veog,:,:));
elseif strcmp(whichchan,'HEOG'); 
    start=2; 
    e=squeeze(D.data(D.channels.heog,:,:));
else
    error('aa:error','\n Oops. Expecting VEOG or HEOG. \n')
end

m=max(abs(e(D.events.start:end,:)),[],1); % get max, abs signal across each trial
s=sort(m);
% (old version looked for steepest smoothed gradient)
% (older version looked for greatest curvature)
% smooth
%[B, A] = butter(1, 2*5/1000,'low');
%f=filtfilt(B,A,s);

if ischar(thresh) % automatically choose threshold maximum signal to allow
%     % OLD METHODS: try to find steepest point or point of maximum curvature
%     % of sorted max values.
%     % get gradient from 5%
%     c=gradient(f);
%     c(1:round(length(c)/20))=0;
%     % choose threshold as point of max gradient
%     [junk threshind]=max(c);
%     thresh=s(threshind);

% NEW METHOD: take threshold (for abs max signal) as 7sd above mean of max(abs(baseline data)).
% (could scale by (whole data length/baseline length)
    bm=max(abs(e(1:D.events.start,:)),[],1);
    thresh=mean(bm) + std(bm)*7;
    %thresh=thresh * sqrt(D.Nsamples/(D.events.start));
    
    maxmethodtext='(auto)';
else
    maxmethodtext='(specified)';
end
r=find(m>thresh);
threshind=find(s>thresh,1);

% reject on basis of mean signal across trial  
% (less senitive to noisy baseline? reject small but stable saccades? early
% saccades penalised more?)  
meansignalacrosstrial=mean(abs(e(D.events.start:end,:)),1);
smean=sort(meansignalacrosstrial);
if ischar(threshmean) % automatically choose threshold  
% (less senitive to noisy baseline? reject small but stable saccades? early saccades penalised more?)    
    meansignalacrossbaseline=mean(abs(e(1:D.events.start,:)),1);
    threshmean=mean(meansignalacrossbaseline)+5*std(meansignalacrossbaseline);
%     threshmeansignal=mean(meansignalacrossbaseline)+7*std(meansignalacrossbaseline);
%     maxsignalacrosstrialsabovethresh=m(meansignalacrosstrial>threshmeansignal);
%     thresh=min(maxsignalacrosstrialsabovethresh);
%     sum(meansignalacrosstrial>threshmeansignal);
%     sum(m>thresh);
    meanmethodtext='(auto)';
else
    meanmethodtext='(specified)';
end
badlist=find(meansignalacrosstrial>threshmean);
threshmeanind=find(smean>threshmean,1);

%% plot diagnostic
subplot(4,2,start+1)
plot(s); hold on; plot(smean,'c');
xlabel('Sorted trials'); ylabel('absolute signal (uV)');axis tight;
%hold on; plot(f,'c'); legend({'Sorted raw data','Low-pass
%filtered'},'Location','Best')
try plot([threshind threshind],ylim,'k--'); catch end
plot(xlim,[thresh thresh],'k--')
try plot([threshmeanind threshmeanind],ylim,'r:'); catch end
plot(xlim,[threshmean threshmean],'r:')
legend({'max across trial','mean across trial'},'location','best');

subplot(4,2,start+2)
hold off; l1=plot(e,'b:'); hold on; axis tight; ylabel('microvolts'); xlabel('Samples');
title({sprintf('%s:\tThreshold_max=%s %s;  %g trials rejected (%.0f%%), leaving %g',whichchan,num2str(thresh),maxmethodtext,length(r),length(r)/length(m)*100,length(m)-length(r)), ...
      sprintf('    \tThreshold_mean=%s %s;  %g trials rejected (%.0f%%), leaving %g',num2str(threshmean),meanmethodtext,length(badlist),length(badlist)/length(m)*100,length(m)-length(badlist))});
l2=plot(e(:,r),'k--');
l3=plot(e(:,badlist),'r:');

try legend([l1(1) l2(1) l3(1)],{'Ok','Rejected by max','Rejected by mean'},'location','best')
catch
    try legend([l1(1) l3(1)],{'Ok','Rejected by mean'},'location','best')
    catch
        try legend([l1(1) l2(1)],{'Ok','Rejected by max'},'location','best')
        catch
        end
    end
end
            
%hold off; plot(c1); axis tight % plot(gradient(gradient(s1))./((1+gradient(s1).^2).^(3/2))); axis tight
colormap jet

return