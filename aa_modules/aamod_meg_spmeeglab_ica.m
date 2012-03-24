% AA module  - intialise a subject by making directory
% Calculate independent components
% Identify and project out blink (and pulse?) components
% Adapted from Jason's _example_pathway.m
% Danny Mitchell 04/04/08

function [aap,resp]=aamod_meg_spmeeglab_ica(aap,task,subjnum,sessnum)

resp='';

switch task
case 'domain'
    resp='session';   % this module needs to be run once per subject

case 'description'
    resp='SPM EEGLAB ICA';
    
case 'summary'
    resp='SPM EEGLAB ICA';
    
case 'checkrequirements'
 
case 'report'

case 'doit'
    
%%%%%%%%%%%%%%% Task-specific setup
subblock=[subjnum sessnum];

% default settings:
I.InputFilter='^S.*<block>.*raw\.mat$';
I.order=5;
I.Overwrite=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check inputs and collect settings if required
% if ~isempty(findstr(char(source),'specify'))
%     if ~isfield(aap.MEG.TaskSettings.(mfilename),'InputFilter')
%         I.InputFilter=spm_input('Filter files to process:',2,'s',I.InputFilter);
%     end
%     if ~isfield(aap.MEG.TaskSettings.(mfilename),'Overwrite')
%         I.Overwrite=spm_input('Overwrite any existing files? ','+1','y/n',[1,0],2);
%     end
% end

aap=aas_meg_filltaskparams(aap, mfilename, I);

%%%%%%%%%%%%% find files and decide whether to run task
if ischar(subblock{1}); return; end
ReportCurrentJob(aap,mfilename,subblock);
files=aas_meg_findfiles(aap,mfilename,subblock);
if isempty(files{1}); fprintf('Found no data! '); return; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do it
addpath /imaging/jt03/spm_eeglab/eeglab_current/eeglab2007November16_beta/;
addpath /imaging/local/spm_eeglab/

for f=1:length(files);
     
    % Run ICA on continuous data (could also do on epoched data):
    % This could take some time, 20 minutes - ? hours?!
    % If you hit memory problems, adjust S.samppct down.
    [pth nam ext]=fileparts(files{f});
    i_file=fullfile(pth,['i' nam ext]);
    if ~exist(i_file,'file') || aap.MEG.TaskSettings.(mfilename).Overwrite
        D=spm_eeg_ldata(files{f});
        clear S
        S.D=files{f};
        S.samppct=1;
        S.newpfx='i';
        try ncomps=D.HarmonicComponents; catch ncomps=D.Nchannels-2; end
        S.args={'extended',1,'pca',ncomps,'maxsteps',800};
        D=spm_eeglab_runica_dm(S); % (don't duplicate data file)
    else load(i_file);
    end
        
    clear S acts
    %%% Get ICA activations:
    Ai_file=fullfile(pth,['A' D.fname]);
    if ~exist(Ai_file,'file') || aap.MEG.TaskSettings.(mfilename).Overwrite
        S.D=i_file;
        S.samppct=1;  % see note below!
        S.newfname=['A' D.fname];
        fprintf('\nCalculating ICs...');
        acts=spm_eeglab_icaact_dm(S); % trying to save memory
        % if adjust samppct, also need to adjust eogdata accordingly:
        % e.g., eogdata=eogdata(:,1:ceil(D.Nsamples*S.samppct));
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
    
    if ~isfield(D.ica,'class'); D.ica.class=''; end;
    if ~isfield(D.ica.class,'blink') || aap.MEG.TaskSettings.(mfilename).Overwrite
        if ~exist('acts','var'); A=spm_eeg_ldata(Ai_file); acts=A.data(:,:); end
        try vdata=D.data(D.channels.veog,:);
        catch D=spm_eeg_ldata(i_file); vdata=D.data(D.channels.veog,:);         
        end
        %%% find blink component...adapted from meg_icablink
        corrthresh=.5; ratiothresh=3;
        fprintf('\nComputing ICA/VEOG correlations...');
        rvec=[];
        for i=1:size(acts,1)
            rvec=[rvec; corr(acts(i,:)',vdata')];
        end
        % Find largest correlations:
        rsort = sort(abs(rvec),'descend');
        rmax=rsort(1); rmax2=rsort(2); rmax3=rsort(3);
        compmax = find(abs(rvec)==rmax);
        compmax2 = find(abs(rvec)==rmax2);
        compmax3 = find(abs(rvec)==rmax3);
        rmaxratio = abs(rvec(compmax))/(abs(rvec(compmax2)));
        fprintf('SUMMARY: 3 largest correlations with VEOG:');
        fprintf('\n%s%d%s%0.3f','IC ',compmax,' = ',rmax);
        fprintf('\n%s%d%s%0.3f','IC ',compmax2,' = ',rmax2);
        fprintf('\n%s%d%s%0.3f','IC ',compmax3,' = ',rmax3);
        fprintf('\n%s%0.3f','--> Ratio 1st:2nd = ',rmaxratio);
        %%% Make blink decision:
        warn=[];
        if rmax>corrthresh
            fprintf('\n%s%s%s','** Corr larger than threshold (',num2str(corrthresh),') **');
            if rmaxratio>ratiothresh
                fprintf('\n%s%s%s','** Ratio larger than threshold (',num2str(ratiothresh),') **');
                fprintf('\n%s%s%s','**** Accepting Component ',num2str(compmax),' as blink component ****');
                blinkcomp=compmax;
            else
                fprintf('\n%s%s%s','!! Ratio smaller than threshold (',num2str(ratiothresh),') !!');
                fprintf('\n%s%s%s','**** Accepting Component ',num2str(compmax),' as blink component ****');
                blinkcomp=compmax;
                warn={'another blink component may exist'};
            end
        else
            fprintf('\n%s%s%s','!! Corr smaller than threshold (',num2str(corrthresh),') !!');
            fprintf('\n%s%s%s','**** Accepting Component ',num2str(compmax),' as blink component ****');
            blinkcomp=compmax;
            warn={'vcorr is low!'};
        end
        if blinkcomp==0; blinkcomp=[]; end
        %%% Store blink component number:
        D.ica.class.blink=blinkcomp;
        D.ica.class.vcorr=rvec(blinkcomp);
        if isfield(D.ica.class,'warnings'); D.ica.class.warnings=[D.ica.class.warnings;warn];
        else D.ica.class.warnings=warn;
        end
        save(fullfile(pth,D.fname),'D');
    end

    try
        if ~isfield(D.ica.class,'pulse') || aap.MEG.TaskSettings.(mfilename).Overwrite
            if ~exist('acts','var'); A=spm_eeg_ldata(Ai_file); acts=A.data(:,:); end
            % find pulse component?
            fprintf('\nEstimating pulse component...');
            apr=zeros(1,size(acts,1));
            hp=spectrum.periodogram('hamming');
            hpopts=psdopts(hp); % hpopts=psdopts(hp,acts(1,:));
            set(hpopts,'Fs',D.Radc,'CenterDC',true,'spectrumtype','onesided');
            pulsecomp=[];
            for i=1:size(acts,1) 
                if i==D.ica.class.blink; continue; end
                Hpsd=psd(hp,acts(i,:),hpopts);
                % psd(hp,log10(Hpsd.data),hpopts);
                apr(i)=(avgpower(Hpsd,[1 1.5])+avgpower(Hpsd,[2.1 2.8])) / ...
                    (avgpower(Hpsd,[1.6 2])+avgpower(Hpsd,[2.8 3.2]));
                 % avgpower(Hpsd,[3.2 3.9]) high?     avgpower(Hpsd,[3.9 4.5]) low?
            end
            [junk pulsecomp]=max(apr);
            fprintf('as #%g ',pulsecomp)
            try fprintf(' (old pulsecomp was #%g)',D.ica.class.pulse); 
                if D.ica.class.pulse~=pulsecomp; figure(3); bar(arp); keyboard; end
            catch
            end
            D.ica.class.pulse=pulsecomp;
            save(fullfile(pth,D.fname),'D');
        end
    catch
        fprintf('failed'); D.ica.class.pulse=[];
        % once spectrum.periodogram('hamming') failed; not sure why;
        % doesn't matter as not projecting out pulse anyway
    end

    %%% Plot topographies and timecourses of selected components, with VEOG:
    [pth na]=fileparts(Ai_file);
    fulloutfile=fullfile(pth,[na '.png']);
    if ~exist(fulloutfile,'file') || aap.MEG.TaskSettings.(mfilename).Overwrite
        if ~exist('acts','var'); A=spm_eeg_ldata(Ai_file); acts=A.data(:,:); end
        myics=union(1:10,[D.ica.class.blink,D.ica.class.pulse]);% Choose ICs you want to plot
        timewin=4*D.Radc; % 4s in samples
        % get data and peak from VEOG timecourse (after 1st event)
        first=D.events.time(1); 
        last=min([D.events.time(end) length(acts)-timewin/2]);
        [junk ind]=max(D.data(D.channels.veog,first:last),[],2);
        ind=ind+first;
        % setup figure
        set(0,'units','pixels'); % this was changed to normalised
        warning off all;
        Fh=spm_figure('GetWin','Graphics'); spm_figure('Clear',Fh);
        warning on all
        pp=get(Fh,'paperposition');
        set(Fh,'paperposition',[pp(1)/2 pp(2)/2 pp(3)+pp(1) pp(4)+pp(2)]);
        colormap jet
        % Plot one at a time:
        spA={}; p=zeros(length(D.channels.eeg),size(acts,1)); % assume all channels
        for i=1:length(myics) % could go from -1 and include EOG at top...
            ic=myics(i);

            sp=subplot('position',[0.03 1-i*(1/length(myics)) 0.17 1/length(myics)]);
            set(sp,'ytick',[],'xtick',[],'xcolor','w','linewidth',2);
            
            if i>0
                % Compute projection:
                clear S proj
                S.D=i_file;
                S.ic2proj=ic;
                S.samppct=.3;  % Don't need much data
                proj=spm_eeglab_icaproj(S);
                p(:,i)=mean(proj,2);
                % plot topography over magnetometers/eeg
                meg_topo_dm(p(:,i),D,'mags',sp)
                ylabel(sprintf('IC %g',ic)); box on
            elseif i==0; xlabel('Magnetometers')
            end
            axis on tight;

            if isempty(regexp(files{f},'-eeg','ONCE'));
                % plot topography over gradiometers
                sp=subplot('position',[0.23 1-i*(1/length(myics)) 0.17 1/length(myics)], ...
                    'ytick',[],'xtick',[],'xcolor','w','linewidth',2);
                if i>0
                    meg_topo_dm(p(:,i),D,'grads',sp); box on
                    if i==ceil(length(myics)/2); ylabel('/\  Magnetometers    |    Gradiometers  \/');end
                elseif i==-1;  ylabel('VEOG')
                elseif i==0; ylabel('HEOG'); xlabel('Magnetometers')
                end
                axis on tight;
            end
        
            % plot 4s period of component timecourse, around peak of VEOG
            spA{i}=subplot('position',[0.4 1-i*(1/length(myics)) 0.6 1/length(myics)]); 
            plot(acts(ic,ind-timewin/2:ind+timewin/2)); axis off tight; 
            % add any events?
            es=find(D.events.time>(ind-timewin/2) & D.events.time<(ind+timewin/2));
            hold on
            plot(repmat(D.events.time(es)-(ind-timewin/2),[2 1]),repmat(ylim',[1,length(es)]),'k:') 
            set(spA{i},'xtick',D.events.time(es)-(ind-timewin/2)); 
            set(spA{i},'xticklabel',D.events.code(es),'xaxislocation','top');
            set(spA{i},'ytick',[],'ycolor',[1 1 1]); box off; axis on
            drawnow
            if ic==D.ica.class.blink % Overlay VEOG timecourse on blink component
                axes('position',get(spA{i},'Position'));
                plot(D.data(D.channels.veog,ind-timewin/2:ind+timewin/2),'r'); axis off tight;
                legend('VEOG signal')
            elseif ic==D.ica.class.pulse % mark in pink
                set(get(gca,'children'),'color','m')
                legend('Pulse component? (not removed)')
            end
        end
        % improve position
        set(Fh,'menuBar','none')
        pos=get(Fh,'position');
        set(Fh,'position',[pos(1) 30 pos(3:4)]);
        % convert topographies to images (gets focus)
        fram=getframe(Fh,[0 0 pos(3)*0.4 pos(4)]); % this gives window focus!
        subplot('position',[0 0 0.4 1]);
        image(fram.cdata); axis off image
        spm_print(strrep(fulloutfile,'png','ps')); % does something to avoid segmentation fault?
        print('-dpng','-r82', fulloutfile); % low res, but filts in web browser
    end

    %%% and project out...
    % perhaps shouldn't try to project out pulse? 
    % See users.fmrib.ox.ac.uk/~rami/fmribplugin
    clear A* acts S
    S.D=i_file;
    S.newfname=['p' D.fname];
    if ~exist(fullfile(pth,S.newfname),'file') || aap.MEG.TaskSettings.(mfilename).Overwrite
        S.ic2rem=[D.ica.class.blink]; %D.ica.class.pulse
        fprintf('\nProjecting out blink component...please wait...\n')
        clear D
        spm_eeglab_icaproj(S);
    end

    % NOTE, you can create a dataset consisting of a single IC
    % component of interest projected to the sensors using:
    % ic2proj=7  % for example
    % clear S
    % S.D=oldfname;
    % S.ic2proj=ic2proj;
    % S.newfname=['proj_' sprintf('%d_',ic2proj) D.fname];
    % [proj D] = spm_eeglab_icaproj(S);
end

otherwise
    aas_log(aap,1,sprintf('Unknown task %s',task));
end;
