function [aap resp]=aamod_emeg_ica_classify(varargin)
% Identify blink (and pulse? and heog?) components
% Adapted from Jason's _example_pathway.m
% Print top 10 components and putative blink and pulse components
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

for f=1:length(files);
clear p

    D=spm_eeg_ldata(files{f}); % ica-ed data
    [pth nam ext]=fileparts(files{f});
    Ai_file=fullfile(pth,['A' nam ext]);
    A=spm_eeg_ldata(Ai_file); % ica activations
    
    if ~isfield(D.ica,'class'); D.ica.class=''; end;

        if length(size(D.data))==2
            vdata=D.data(D.channels.veog,:);
            hdata=D.data(D.channels.heog,:);
            acts=A.data(:,:);
        else
            % assuming we want to exlude rejected trials
            epindex=1:D.Nevents;
            vdata=squeeze(D.data(D.channels.veog,:,epindex(~D.events.reject)));
            hdata=squeeze(D.data(D.channels.heog,:,epindex(~D.events.reject)));
            vdata=reshape(vdata,1,size(vdata,1)*size(vdata,2));
            hdata=reshape(hdata,1,size(hdata,1)*size(hdata,2));
            acts=A.data(:,:,:);
            acts=reshape(acts,size(acts,1),size(acts,2)*size(acts,3));
        end
    
%% find possible eye movement components...adapted from meg_icablink
    if ~isfield(D.ica.class,'blink') || ~isfield(D.ica.class,'heog') || ...
            settings.Overwrite
        
        %corrthresh=.5; ratiothresh=3;
        fprintf('\nComputing correlations with VEOG/HEOG...');
        Vrvec=[]; Hrvec=[];
        for i=1:size(acts,1)
            Vrvec=[Vrvec; corr(acts(i,:)',vdata')];
            Hrvec=[Hrvec; corr(acts(i,:)',hdata')];
        end
        % Find above-threshold correlations and
        % store blink and saccade component numbers:
        D.ica.class.blink=find(abs(Vrvec)>settings.VeogCorrThresh);
        D.ica.class.vcorr=Vrvec; % store all correlation values
        
        D.ica.class.heog=find(abs(Hrvec)>settings.HeogCorrThresh);        
        D.ica.class.hcorr=Hrvec;        
    end
    
%% find possible pulse component
    if ~isfield(D.ica.class,'pulse') || settings.Overwrite
        if ~exist('p','var'); p=getspatialpatterns(D,acts,files{f},settings.Overwrite); end;
        fprintf('\nEstimating pulse components');
        apr=zeros(size(acts,1),2);
        patdif=zeros(1,size(acts,1));
        %% get example pulse spike, assuming samples at 250Hz and pulse
        %% rate < 2Hz, and pattern over sensors
        load /imaging/dm01/MEG/aaMEG/pulsespike_250samples;
        load /imaging/dm01/MEG/aaMEG/pulsepattern_mean;
        temp=[];
        for i=1:size(acts,1)
            if ~isempty(intersect(i,union(D.ica.class.blink,D.ica.class.heog)));
                apr(i,:)=[0 0];
            else
                % temporal: convolve with spike detector
                z=zscore(acts(i,:)); % make it fairer?
                mysigns=[1 -1];
                for ms=1:2
                spikes=conv(z,mysigns(ms)*pulsespike); % this was absed before;
                spikes=sort(spikes,'descend');
                %temp=[temp; spikes];
                apr(i,ms)=sum(spikes(1000:10000)); 
                % range selected from diagnostic plot below, for 15min of
                % data; might vary by data length, so proportion might be
                % better...
                end
                
                % spatial: 
                if strcmp(D.modality,'MEG')
                    patdif(i)=sum((zscore(abs(p(:,i)))-abs(pulsepattern)).^2);
                else patdif(i)=1; % pulsepattern is defined for 306 meg channels, not eeg electrodes
                end
            end
            fprintf('.')
        end
        D.ica.class.pulsescores=range(zscore(apr),2)'-zscore(patdif).*double(all(apr')>0);       
        [junk maxind]=max(D.ica.class.pulsescores);
        D.ica.class.pulse=intersect(find(D.ica.class.pulsescores>settings.PulseThresh),maxind);   
        
%         tits={'Pulseness';'PulseRat'};
%         vals={toppulsescore;toppulsescore/nextpulsescore};
%         defs={'Highest estimate of how much IC looks like pulse, based on temporal and spatial patterns';...
%             'Ratio of highest and second highest pulseness scores'};
%         aas_emeg_savestats(D,tits,vals,defs);
% 
%         %% temp: for constructing average pulse pattern
%         try
%             temppulsefilename=fullfile('/imaging/dm01/MEG/aaMEG',['pulsepat_' regexprep(D.path,'.*/(meg.......)/','$1') D.fname]);
%             temppulsepat=p(:,pulsecomp);
%             save(temppulsefilename,'temppulsepat')
%         catch
%         end
%         %%       
        %%% diagnostic plot
        % plot(temp');axis tight; hold on; set(gca,'xScale','log'); 
        % plot(temp(pulsecomp,:),'ko');
        
        fprintf('as #%s ',mat2str(D.ica.class.pulse))
    end
    save(fullfile(pth,D.fname),'D');

%% Prepare to plot topographies and timecourses of selected components, with EOG:
    [pth na]=fileparts(Ai_file);
    fulloutfile=fullfile(pth,'figures',[na '.png']);
    if ~exist(fulloutfile,'file') || settings.Overwrite
        if ~exist('acts','var'); A=spm_eeg_ldata(Ai_file); acts=A.data(:,:); end
         if ~exist('p','var'); p=getspatialpatterns(D,acts,files{f}); end;
         % Choose ICs to plot
         try
            %myics=union(1:9,[D.ica.class.blink', D.ica.class.pulse',D.ica.class.pulse2', D.ica.class.heog']);
            marked=unique([D.ica.class.blink(:); D.ica.class.pulse(:); D.ica.class.heog(:)]);
            rest=setxor(1:size(acts,1),marked);
            myics=sort([marked' rest(1:(8-length(marked)))]);
         catch
             debugnow
         end
        timewin=4*D.Radc; % 4s in samples
        % get data and peak from EOG timecourses (after 1st event)
%         first=D.events.time(1);
%         last=min([D.events.time(end) length(acts)-timewin/2]);
%         [junk ind]=max(D.data(D.channels.veog,first:last),[],2);
%         ind=ind+first;
%         [junk Hind]=max(D.data(D.channels.heog,first:last),[],2);
%         Hind=Hind+first;
        set(0,'units','pixels'); % this was changed to normalised
        warning off all; Fh=spm_figure('GetWin','Graphics'); spm_figure('Clear',Fh); warning on all
        pp=get(Fh,'paperposition');
        set(Fh,'paperposition',[pp(1)/2 pp(2)/2 pp(3)+pp(1) pp(4)+pp(2)],'renderer','painters'); 
        set(Fh,'menuBar','none')
        
%% plot 4s period of component timecourse, around peak of component, or EOG,
% from middle of timecourse

        [junk Vind]=max(vdata(timewin:end-timewin),[],2);
        [junk Hind]=max(hdata(timewin:end-timewin),[],2);
        spA={}; H1=[]; H2=[]; H3=[];    
        for i=1:length(myics)
            ic=myics(i);
            if ismember(ic,D.ica.class.blink); ind=Vind; 
            elseif ismember(ic,D.ica.class.heog); ind=Hind; 
            else [junk ind]=max(acts(ic,timewin:end-timewin),[],2);
            end
            ind=ind+timewin;
            spA{i}=subplot('position',[0.4 .98-i*.1 0.6 .1]);
            try
                plot(acts(ic,ind-timewin/2:ind+timewin/2)); axis off tight;
            catch
                fprintf('\nPlease debug %s\n',mfilename); keyboard
            end
            % add any events?
            if length(size(D.data))==2
                eventtimes=D.events.time;
                eventcodes=D.events.code;
            else
                eventtimes=find(mod(1:size(acts,2),D.Nsamples)==D.events.start);
                eventcodes=A.events.code; % (rejected trials may have been removed)
            end
            es=find(eventtimes>(ind-timewin/2) & eventtimes<(ind+timewin/2));
            hold on
            plot(repmat(eventtimes(es)-(ind-timewin/2),[2 1]),repmat(ylim',[1,length(es)]),'k:')
            set(spA{i},'xtick',eventtimes(es)-(ind-timewin/2));
            set(spA{i},'xticklabel',eventcodes(es),'xaxislocation','top');
            set(spA{i},'ytick',[],'ycolor',[1 1 1]); box off; axis on
            drawnow
            if ismember(ic,D.ica.class.pulse) % mark putative pulse component in pink
                H3=get(gca,'children'); H3=H3(end); set(H3,'color','m');
                text(0,max(ylim)/2,sprintf('pulsescore=%.2g',D.ica.class.pulsescores(ic)),'fontweight','bold')
            end
            if ismember(ic,D.ica.class.blink) % Overlay VEOG timecourse on putative blink component
                axes('position',get(spA{i},'Position'));
                H1=plot(vdata(ind-timewin/2:ind+timewin/2),'r'); axis off tight;
                text(0,max(ylim)/2,sprintf('r=%.2g',D.ica.class.vcorr(ic)),'fontweight','bold')
            end
            if ismember(ic,D.ica.class.heog) % Overlay HEOG timecourse
                axes('position',get(spA{i},'Position'));
                H2=plot(hdata(ind-timewin/2:ind+timewin/2),'r:'); axis off tight;
                text(0,max(ylim)/2,sprintf('r=%.2g',D.ica.class.hcorr(ic)),'fontweight','bold')
            end
        end
        % add legend
        try
            labs={'VEOG','HEOG','Pulse?'};
            legend([H1 H2 H3(end)],labs(~cellfun('isempty',{H1 H2 H3})), ...
                'orientation','horizontal','fontsize',8);
        catch 
        end
        
%% Plot projections to sensors:
        fprintf('\nPlotting topographies')
        load /imaging/dm01/MEG/aaMEG/redblue2_128
        colormap(map)
        for r=1:length(myics) % could go from -1 and include EOG at top...
            ic=myics(r);

            sp=subplot('position',[0.03 .98-r*.1 0.17 .1]);
            
            meg_topo_dm(p(:,ic),D,'mags',sp,0); % plot topography over magnetometers/eeg
            % (final flag specifies image bitmap rather than vector graphics)
            set(sp,'ytick',[],'xtick',[],'xcolor','w','linewidth',2);
            ylabel(sprintf('IC %g',ic),'fontweight','bold');
            box on
            axis on tight;

            if isempty(regexp(files{f},'-eeg','ONCE'));
                sp=subplot('position',[0.23 .98-r*.1 0.17 .1]);
                    meg_topo_dm(p(:,ic),D,'grads',sp,0); box on % plot topography over gradiometers
                    if i==floor(length(myics)/2); 
                        ylabel('/\  Magnetometers    |    Gradiometers  \/', 'fontsize', 12);
                    end
                set(sp,'ytick',[],'xtick',[],'xcolor','w','linewidth',2);
                axis on tight;
            end
            fprintf('.'); drawnow
        end
        
        %%% Plot Z-score histogram (nice idea from Jason):
        xlabs={'correlation with VEOG','correlation with HEOG','pulseness'};
        for a=1:3
            subplot('position',[(a-1)/3+0.04 .06 0.29 0.1]);
            switch a
                case 1; his=zscore(D.ica.class.vcorr); marked=D.ica.class.blink;
                case 2; his=zscore(D.ica.class.hcorr); marked=D.ica.class.heog;
                case 3; his=D.ica.class.pulsescores; marked=D.ica.class.pulse;
            end
            hist(his,40)
            set(gca,'yTick',[]); box off
            xlabel(sprintf('Z-score of %s',xlabs{a}));
            if a==1; ylabel('Frequency'); end

            % Label high outliers:
            for m=1:length(marked)
                text(his(marked(m)),4,sprintf('IC %d',marked(m)),'horizontalalign','center','fontweight','bold');
            end
        end
    
        spm_print(strrep(fulloutfile,'png','ps')); % does something to avoid segmentation fault?
        print('-dpng','-r82', fulloutfile); % low res, but fits in web browser
    end

end % next file

return

%% get spatial patterns
function p=getspatialpatterns(D,acts,i_file,overwrite)

fname=fullfile(D.path,['ICs_' D.fname]);
if exist(fname,'file') && ~overwrite; 
    fprintf('\nLoading spatial projections') 
    load(fname);
else
    p=zeros(length(D.channels.eeg),size(acts,1)); % assume all channels
    fprintf('\nComputing and saving spatial projections')
    clear S proj
    S.D=i_file;
    S.samppct=.2;  % Don't need much data
    for i=1:size(acts,1)
        S.ic2proj=i;
        proj=spm_eeglab_icaproj(S);
        p(:,i)=mean(proj,2);
        fprintf('.')
    end
    save(fname,'p');
end

return