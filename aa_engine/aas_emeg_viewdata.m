function [varargout]=aas_emeg_viewdata(varargin)
% Danny Mitchell 2008

%% get inputs
if length(varargin)>0
    D=varargin{1};
else
    clear
    [D ok]=spm_select(1,'mat','Please select a data file','',pwd,'^.*mat$');
    if ~ok; return;end
end

global e map continuous epoched erf tc start clim chanstuff Cpos meg H magname gradname
varargout={};

%% load data
fprintf('\nLoading data...')
if ischar(D)
    D=spm_eeg_ldata(char(D));
elseif ~isfield(D,'data') || isempty(D.data)
    D=spm_eeg_ldata(fullfile(D.path,D.fname));
end
D.data=D.data(:,:,:);

%% determine data type
continuous=false; epoched=false; erf=false;
if ndims(D.data)==2; continuous=true;
elseif D.Nevents==length(D.events.types); erf=true;
else epoched=true;
end

if strcmp(D.modality,'MEG'); meg=true;
else meg=false;
end

%% create erfs
if epoched
    epochs=D.data;
    D.data=zeros(size(D.data,1),size(D.data,2),length(D.events.types));
    for t=1:length(D.events.types)
        toave=find(D.events.code==D.events.types(t) & ~D.events.reject);
        D.data(:,:,t)=mean(epochs(:,:,toave),3);
        % EOG 95% confidence intervals (uncorrected)
        [h p ci]=ttest(squeeze(epochs(D.channels.veog,:,toave))',0,0.05,'both');
        D.Vci(:,:,t)=ci;
        [h p ci]=ttest(squeeze(epochs(D.channels.heog,:,toave))',0,0.05,'both');
        D.Hci(:,:,t)=ci;
        fprintf('.')
    end

end
fprintf('Done')

%% define some inital parameters
Sm=Inf; Sg=Inf;% will contain handles to selected channel markers
magname=D.channels.name(1);
if meg; gradname={'MEG0111'};
else gradname={};
end
if continuous; width=4; % seconds per window
    step=fix(width*D.Radc); % seconds to samples per window
else step=D.Nsamples;
end
tc=fix(step/2)+1; % sample of topoplot; at centre of window for continuous data
e=1; % event type (not relevant for continuous)

%% initialise figure
fprintf('\nInitialising figure')
F=figure(999); clf
set(0,'units','pixels');
set(F,'paperunits','normalized','papersize',[1 1],'paperType','<custom>', ...
    'paperPosition',[0 0 1 1],'paperOrientation','landscape')
set(F,'position',[6 32 1200 840],'menubar','none');
set(F,'name',fullfile(D.path,D.fname));
addpath /imaging/dm01/MEG/aaMEG
load /imaging/dm01/MEG/aaMEG/redblue2; % colormap in map

%% Channel stuff for topographies:
try load(D.channels.ctf);
catch % might be EEG template file?
    load(fullfile('/imaging/local/spm/spm5/EEGtemplates',D.channels.ctf));
end
fprintf('.')

% Create channel matrix 'ch':
%  locations(1:102) X sensor type (1,2,3)
names=D.channels.name;
chanstuff.ch=[];
if length(D.channels.eeg)==306 % djm
    for i=1:102;
        chanstuff.ch=[chanstuff.ch;find(strncmpi(names,names{i}(1:6),6))'];
    end
else chanstuff.ch=(1:length(D.channels.eeg))';
end

chanstuff.Cpos = Cpos(:, D.channels.order(1:length(chanstuff.ch))); % djm (was 102)
chanstuff.x = min(chanstuff.Cpos(1,:)):0.01:max(chanstuff.Cpos(1,:)); % halved default res to hopefully speed things up
chanstuff.y = min(chanstuff.Cpos(2,:)):0.01:max(chanstuff.Cpos(2,:));
[chanstuff.x1,chanstuff.y1] = meshgrid(chanstuff.x,chanstuff.y);
chanstuff.xp = chanstuff.Cpos(1,:)';
chanstuff.yp = chanstuff.Cpos(2,:)';
chanstuff.Cnames=Cnames;

[LHS RHS]=getLRpairs2(D,'sideonly_noplot');

%% initialise plots - this should hopefully save time later
A=[];
H.gfp=subplot('position',           [0.05 0.05 0.7 0.18]);
FL=ylabel('Global field power'); A=lab2ann(gca,FL,A);
if continuous
    FL=xlabel('Time (s). Recentre at click, or click left margin to select time...');A=lab2ann(gca,FL,A);
else
    FL=xlabel('Time (ms). Reset topography at click, or click left margin to select time...');A=lab2ann(gca,FL,A);
end
H.selection=subplot('position',     [0.05 0.23 0.7 0.18]);
if meg
    H.butterflygrads=subplot('position',[0.05 0.41 0.7 0.18]);
    FL=ylabel('Absolute flux gradient (a.u.)');A=lab2ann(gca,FL,A);
end
H.butterflymags=subplot('position', [0.05 0.59 0.7 0.18]);
if meg; FL=ylabel('Magnetic flux (fT)');
else FL=ylabel('Potential (uV)');
end
A=lab2ann(gca,FL,A);
H.eog=subplot('position',           [0.05 0.77 0.7 0.18]);
FL=ylabel('EOG (uV)');A=lab2ann(gca,FL,A); hold on

plot(1:2,'linewidth',2,'color',[0.5 0 1]);
plot(1:2,'linewidth',2,'color',[0.5 0.5 0]);
plot(1:2,'linewidth',2,'color',[0 0 1]);
plot(1:2,'linewidth',2,'color',[1 0 0]);
plot(1:2,'linewidth',2,'color',[0 0 0.7]);
plot(1:2,'linewidth',2,'color',[0 0 0]);
plot(1:2,'linewidth',2,'color',[0 0.5 0]);
if meg; leglabs={'VEOG','HEOG','Left sensors','Right sensors','Selected magnetometers','Selected gradiometer pairs','Global field power'};
else leglabs={'VEOG','HEOG','Left sensors','Right sensors','Selected electrodes','Selected - contralateral','Global field power'};
end
L=legend(leglabs, 'location',[0.8 0.81 0.12 0.14]); copyobj(L,F);
fprintf('.')

H.magtopo=subplot('position',  [0.76 0.50 0.23 0.3]);
set(H.magtopo,'yaxislocation','right');
if meg; FL=ylabel('Magnetic flux');
else FL=ylabel('Scalp potential');
end
A=lab2ann(gca,FL,A);
topoplot(D.data(:,1,e),'mags');cla
tS=plot(chanstuff.xp(:)*size(chanstuff.x1,2),(chanstuff.yp(:))*size(chanstuff.x1,1),'w.');
H.Sm=copyobj(H.magtopo,F); axes(H.Sm); axis off; delete(tS)

if meg
    H.gradtopo=subplot('position', [0.76 0.20 0.23 0.3]);
    set(H.gradtopo,'yaxislocation','right');
    FL=ylabel('Absolute magnetic flux gradient');A=lab2ann(gca,FL,A);
    FL=xlabel({'Click to select nearest channel','or click and drag to average multiple channels'}); lab2ann(gca,FL,A);
    topoplot(D.data(:,1,e),'grads');cla
    tS=plot(chanstuff.xp(:)*size(chanstuff.x1,2),(chanstuff.yp(:))*size(chanstuff.x1,1),'w.');
    H.Sg=copyobj(H.gradtopo,F); axes(H.Sg); axis off; delete(tS)
else H.Sg=-1; H.gradtopo=-1;
end

if ~continuous
    H.previous=uicontrol('style','pushbutton','units','normalized','position',[0.05 0.002 0.086 0.025], ...
        'string','Previous event');
    H.next=uicontrol('style','pushbutton','units','normalized','position',[0.664 0.002 0.086 0.025], ...
        'string','Next event');
else H.previous=[];H.next=[];
end

H.print=uicontrol('style','pushbutton','units','normalized','position',[0.77 0.102 0.22 0.025], ...
    'string','Print figure');
H.copywave=uicontrol('style','pushbutton','units','normalized','position',[0.77 0.077 0.22 0.025], ...
    'string','Collect selected timecourses');
H.copytopo=uicontrol('style','pushbutton','units','normalized','position',[0.77 0.052 0.22 0.025], ...
    'string','Collect current topographies');
H.return=uicontrol('style','pushbutton','units','normalized','position',[0.77 0.002 0.22 0.025], ...
    'string','Exit and return data structure');

fprintf('Done\n')

%% do it
[Sm Sg]=doplots(D,step,LHS,RHS,Sm,Sg);
while 1==1

    try axes(H.previous); axes(H.next); catch end

    try [x y]=ginput(1);
    catch fprintf('Exited.\n'); return;
    end

    if gco==H.print
        [pth nam]=fileparts(D.fname);
        fout=inputdlg('Enter file name:','Append figure to postscript file',[1 99],{fullfile(D.path,[nam '_fig.ps'])});
        if ~isempty(fout)
            try
                M=msgbox(sprintf('Printing to: \n%s',fout{1}),'Please wait...');
                print(999,'-dpsc2',fout{1},'-r200','-append','-noui');
                delete(M)
            catch
                keyboard
                delete(M)
                M=msgbox(sprintf('Failed to print to %s',fout{1}),'Warning','warn','modal');
            end
        end
    elseif gco==H.previous; e=e-1; [Sm Sg]=changeevent(D,Sm,Sg);
    elseif gco==H.next; e=e+1; [Sm Sg]=changeevent(D,Sm,Sg);
    elseif gco==H.copytopo; copytopo;
    elseif gco==H.copywave; copywave(D);
    elseif gco==H.return; varargout{1}=D; return
    elseif x<min([0 start]); inp=[]; % left margin to specify time
        while isempty(inp)
            inp=inputdlg('Jump to time (s):');
            if isempty(inp); break; end
            try
                switch inp{1};
                    case {'blink','veog'};
                        [junk skip]=max(D.data(D.channels.veog,tc:end,e));
                        tc=tc+skip;
                    case {'saccade','heog'};
                        [junk skip]=max(abs(diff(D.data(D.channels.heog,tc:end,e))));
                        tc=tc+skip;
                    case {'max','peak'};
                        [junk skip]=max(mean(abs(D.data(D.channels.eeg,tc:end,e))));
                        tc=tc+skip;
                    otherwise
                        inp=eval(inp{1});
                        tc=fix(inp*D.Radc);
                end
            catch inp=[];
            end
        end
    elseif gca==H.Sm || gca==H.magtopo % select channel
        if exist('Sm','var') && all(ishandle(Sm)); delete(Sm); end
        [selectedm magname]=selectchans(x,y);
    elseif gca==H.Sg || gca==H.gradtopo
        % select channel
        if exist('Sg','var') && all(ishandle(Sg)); delete(Sg); end
        [selectedg gradname]=selectchans(x,y);
    else % select timepoint(s)
        %if gca==H.eog; beep;beep;beep; keyboard; end % FOR DEBUGGING!!!
        %tc=start+fix(x);
        rbbox;                   % return figure units
        point2 = get(gca,'CurrentPoint');    % button up detected
        x2 = point2(1,1);
        tc=start+fix([min([x,x2]) max([x,x2])]);
    end
    if tc<1; tc=1; end
    if tc>D.Nsamples; tc=D.Nsamples; end
    [Sm Sg]=doplots(D,step,LHS,RHS,Sm,Sg);
    drawnow
end

%%
function [Sm Sg]=changeevent(D,Sm,Sg)
global e
if e<1; e=D.events.Ntypes;
elseif e>D.events.Ntypes; e=1;
end
delete(Sm);
if ishandle(Sg); delete(Sg); end
Sm=Inf; Sg=Inf;% S=Inf forces replot of all axes
return

%%
function [selected selectednames]=selectchans(x,y)
global chanstuff
xp=chanstuff.xp*size(chanstuff.x1,2);
yp=chanstuff.yp*size(chanstuff.y1,1);
rbbox;                   % return figure units
point2 = get(gca,'CurrentPoint');    % button up detected
[x2 y2] = deal(point2(1,1),point2(1,2));              % extract x and y
selected=find(xp>min(x,x2) & xp<max(x,x2) & yp>min(y,y2) & yp<max(y,y2));
if isempty(selected) % get closest to click
    dist=sqrt((xp-x).^2+(yp-y).^2);
    [junk selected]=min(dist);
end
selectednames=cellstr(char(chanstuff.Cnames{selected})); % mag/eeg name(s) at this position
return

%% update all the plots
function [Sm Sg]=doplots(D,step,LHS,RHS,Sm,Sg)
global e tc continuous epoched start clim chanstuff meg H magname gradname
ticksteps=16;

start=tc(1)-fix(step/2);
stop=tc(1)+fix(step/2)-1;
if start<1; start=1; stop=step; end
if stop>D.Nsamples; stop=D.Nsamples; start=D.Nsamples-step+1; end
data=D.data(:,start:stop,e);

%% clear time markers
tm=findobj(999,'userdata','timemarker');
delete(tm);

%% get any triggers and set window title
if continuous
    trigs=find(D.events.time>(tc(1)-fix(step/2)) & D.events.time<(tc(1)+fix(step/2)));
    trigtimes=D.events.time(trigs);
    trigtimes=trigtimes-start; % convert to axis space
    trigcodes=D.events.code(trigs);
else
    trigtimes=D.events.start; trigcodes=[];
    % why D.events.start+1 and not D.events.start?? am I out by a sample??
    try enam=D.events.names{e}; catch enam=''; end
    set(999,'name',fullfile(D.path,sprintf('%s - Event %g of %g: %s',D.fname,e,length(D.events.types),enam)));
end

%% eog, gfp and butterfly plots
if ~continuous && ~all(isinf(Sm))
    % plots initialised; no need to replot all timecourses unless continuous
    % but need to update topo marker...
    axes(H.butterflymags); waveplot();
    if meg; axes(H.butterflygrads); waveplot(); end
    axes(H.eog); waveplot();
    axes(H.gfp); waveplot();
else
    if meg; sensperside=48;
        g=sqrt(data([LHS(49:96) RHS(49:96)],:).^2+data([LHS(97:144) RHS(97:144)],:).^2);
    else sensperside=length(LHS);
    end
    if ~continuous
        if meg
            clim.mags=[min(min(data(1:102,:))) max(max(data(1:102,:)))];
            d=sqrt(data(chanstuff.ch(:,2),:).^2 + data(chanstuff.ch(:,3),:).^2);
            clim.grads=[min(d(:)) max(max(d(:)))];
        else
            clim.mags=[min(min(data(D.channels.eeg,:))) max(max(data(D.channels.eeg,:)))];
        end

        axes(H.butterflymags); waveplot(data(LHS(1:sensperside),:)',trigtimes,{'b'},'nomarker');
        axes(H.butterflymags); waveplot(data(RHS(1:sensperside),:)',trigtimes,{'r'},'holdon');
        if meg
            axes(H.butterflygrads); waveplot(g(1:sensperside,:)',trigtimes,{'b'},'nomarker');
            axes(H.butterflygrads); waveplot(g((sensperside+1):end,:)',trigtimes,{'r'},'holdon');
        end
    else
        if meg % to speed things up a bit
            sensperside=20;
            [junk indL]=sort(range(g(1:48,:),2),'descend');
            [junk indR]=sort(range(g(49:end,:),2),'descend');
            axes(H.butterflygrads); waveplot(g(indL(1:sensperside),:)',trigtimes,{'b'},'nomarker');
            axes(H.butterflygrads); waveplot(g(indR(1:sensperside)+48,:)',trigtimes,{'r'},'holdon');

            [junk indL]=sort(range(data(LHS(1:48)),2),'descend');
            [junk indR]=sort(range(data(RHS(1:48)),2),'descend');
        else
            indL=1:length(LHS);
            indR=1:length(RHS);
        end
        axes(H.butterflymags); waveplot(data(LHS(indL(1:sensperside)),:)',trigtimes,{'b'},'nomarker');
        axes(H.butterflymags); waveplot(data(RHS(indR(1:sensperside)),:)',trigtimes,{'r'},'holdon');
    end
    set(H.butterflymags,'xtick',[]);
    if meg; set(H.butterflygrads,'xtick',[]); end

    axes(H.eog);
    waveplot(data(D.channels.veog,:)',[],{'color',[0.5 0 1]},'nomarker');
    waveplot(data(D.channels.heog,:)',trigtimes,{'color',[0.5 0.5 0]},'holdon');
    if epoched
        waveplot(D.Vci(:,:,e)',[],{'color',[0.5 0 1],'linestyle',':'},'extra');
        waveplot(D.Hci(:,:,e)',[],{'color',[0.5 0.5 0],'linestyle',':'},'extra');
        ylims=ylim; r=range(ylims);
        waveplot((D.Vci(1,:,e)>0 | D.Vci(2,:,e)<0)*r+ylims(1)-r/12,[],{'color',[0.5 0 1],'linestyle','none','marker','.'},'extra');
        waveplot((D.Hci(1,:,e)>0 | D.Hci(2,:,e)<0)*r+ylims(1)-2*r/12,[],{'color',[0.5 0.5 0],'linestyle','none','marker','.'},'extra');
    end
    set(H.eog,'xaxisLocation','top','xtick',trigtimes,'xticklabel',trigcodes,'tickLength',[0 0]);

    axes(H.gfp); waveplot(mean(abs(data)), trigtimes,{'color',[0 0.5 0],'linewidth',1});
end

%% plot selected sensors
axes(H.selection);
%selectedm=find(strcmpi(D.channels.name,magname));
selectedm=find(ismember(D.channels.name,magname));
waveplot(mean(data(selectedm,:),1),trigtimes,{'color',[0 0 0.7],'linewidth',1});
if meg
    %selectedg=find(strncmpi(D.channels.name,gradname(1:6),6));
    sg2=find(ismember(D.channels.name,regexprep(gradname,'1$','2')));
    sg3=find(ismember(D.channels.name,regexprep(gradname,'1$','3')));
    waveplot(mean(sqrt(data(sg2,:).^2+data(sg3,:).^2),1),trigtimes,{'color',[0 0 0],'linewidth',1},'holdon');
    if length(magname)==1; ylabel(sprintf('%s & %s',magname{1},gradname{1})); else ylabel(''); end
else
    if length(magname)==1; ylabel(sprintf('Electrode: %s',magname{1})); else ylabel(''); end
    selleft=ismember(selectedm,LHS);
    selright=ismember(selectedm,RHS);
    if any(ismember(selectedm,LHS))
        contra=D.channels.eeg(RHS(ismember(LHS,selectedm(selleft))));
        ipsi=D.channels.eeg(LHS(ismember(LHS,selectedm(selleft))));
    else contra=[]; ipsi=[];
    end
    if any(ismember(selectedm,RHS))
        contra=[contra D.channels.eeg(LHS(ismember(RHS,selectedm(selright))))];
        ipsi=[ipsi D.channels.eeg(RHS(ismember(RHS,selectedm(selright))))];
    end
    if ~isempty(contra)
        waveplot(mean(data(ipsi(:),:)-data(contra(:),:),1),trigtimes,{'color',[0 0 0],'linewidth',1},'holdon');
    else % probably midline electrode
    end
end
set(H.selection,'xtick',[]);

%% finalise axis labels
xlims=xlim;
ticks=xlims(1):step/ticksteps:xlims(2);
if ~continuous; ticks=sort(unique([ticks xlims(2)])); end

axes(H.gfp)
if continuous; labels=(start:step/ticksteps:(stop-1))/D.Radc; % s from start of acquisition
else
    labels=((-D.events.start:step/ticksteps:D.events.stop)+1)/D.Radc*1000; % ms from trigger
    labels=sort(unique([labels D.events.stop/D.Radc*1000]));
end
set(gca,'xtick',ticks,'xticklabel',labels );

%% plot topographies and mark selected sensors
axes(H.magtopo); topoplot(mean(D.data(:,tc(1):tc(end),e),2),'mags');
if continuous; 
    if tc(1)==tc(end); xlabel(sprintf('%g ms',tc(1)/D.Radc*1000));
    else xlabel(sprintf('%g-%g ms',tc(1)/D.Radc*1000,tc(2)/D.Radc*1000));
    end
else
    if tc(1)==tc(end); xlabel(sprintf('%g ms',(tc(1)-D.events.start-1)/D.Radc*1000));
    else xlabel(sprintf('%g-%g ms',(tc(1)-D.events.start-1)/D.Radc*1000,(tc(2)-D.events.start-1)/D.Radc*1000));
    end     
end

axes(H.Sm); % put sensors back on top
if ~all(ishandle(Sm));
    clear Sm
    Sm(1)=plot(chanstuff.xp(selectedm)*size(chanstuff.x1,2),(chanstuff.yp(selectedm))*size(chanstuff.x1,1),'k.');
    Sm(2)=plot(chanstuff.xp(selectedm)*size(chanstuff.x1,2),(chanstuff.yp(selectedm))*size(chanstuff.x1,1),'ow');
end

if meg
    axes(H.gradtopo); topoplot(mean(D.data(:,tc(1):tc(end),e),2),'grads');
if continuous; 
    if tc(1)==tc(end); xlabel(sprintf('%g ms',tc(1)/D.Radc*1000));
    else xlabel(sprintf('%g-%g ms',tc(1)/D.Radc*1000,tc(2)/D.Radc*1000));
    end
else
    if tc(1)==tc(end); xlabel(sprintf('%g ms',(tc(1)-D.events.start-1)/D.Radc*1000));
    else xlabel(sprintf('%g-%g ms',(tc(1)-D.events.start-1)/D.Radc*1000,(tc(2)-D.events.start-1)/D.Radc*1000));
    end     
end
    axes(H.Sg);
    if ~all(ishandle(Sg));
        clear Sg
        sg1=find(ismember(D.channels.name,regexprep(gradname,'2/3','1')));
        Sg(1)=plot(chanstuff.xp(sg1)*size(chanstuff.x1,2),(chanstuff.yp(sg1))*size(chanstuff.x1,1),'k.');
        Sg(2)=plot(chanstuff.xp(sg1)*size(chanstuff.x1,2),(chanstuff.yp(sg1))*size(chanstuff.x1,1),'ow');
        % note sg1 are actually the magnetometers at the positions of the
        % selected gradiometers
    end
end



return

%%
function waveplot(data,trigtimes,style,note)
global tc start

if ~exist('note','var'); note=''; end
if ~strcmp(note,'holdon') && ~strcmp(note,'nomarkerholdon') && ~strcmp(note,'extra'); hold off; end

if ~exist('data','var') % update time marker
    hold on
else
    if exist('style','var'); plot(data,style{:})
    else plot(data)
    end
    hold on
    if ~strcmp(note,'extra');
        axis tight;
        if ~isempty(trigtimes)
            plot(repmat(trigtimes,[2 1]),repmat(ylim',[1,length(trigtimes)]),'k--')
        end
    end
end

if ~strcmp(note,'nomarker') && ~strcmp(note,'nomarkerholdon') && ~strcmp(note,'extra')
    for t=1:length(unique(tc))
        tm=plot(repmat(tc(t)-start,[2 1]),ylim,'color',[0.9 0.9 0.9]);
        set(tm,'userdata','timemarker');
    end
end

return

%%
function topoplot(td,chtype)
cla
global map clim continuous chanstuff

switch chtype
    case {'mags','eeg'}
        d=td(chanstuff.ch(:,1));
    case {'grads'}
        % Get grad vector length:
        d2=td(chanstuff.ch(:,2));
        d3=td(chanstuff.ch(:,3));
        d=sqrt(d2.^2 + d3.^2);
end

z = griddata(chanstuff.xp, chanstuff.yp, d, chanstuff.x1, chanstuff.y1, 'cubic');

if ~continuous && ~isempty(clim) && isstruct(clim) && isfield(clim,chtype)
    z=z-clim.(chtype)(1);
    z=z/diff(clim.(chtype))*64;
else % use full colormap
    z=z-min(z(:));
    z=z/max(z(:))*64;
end

switch chtype
    case {'mags','eeg'}
        try subimage(z,colormap(map))
        catch subimage(z,colormap(jet))
        end
    case {'grads'}
        subimage(z,colormap(hot))
end

axis xy image
set(gca,'xtick',[],'ytick',[])
hold on
return

%%
function copytopo
global H clim map meg
% colourbars need some work!

figure(991)
% get current plots
M=findobj(991,'type','axes','tag','mags');
G=findobj(991,'type','axes','tag','grads');
C=findobj(991,'tag','colorbar');
delete(C);

% reposition
if ~isempty(clim) && isstruct(clim) && isfield(clim,'grads')
    w=0.9;% leave space for colorbar
else w=1;
end

%r=ceil(sqrt(length(A)+1));
c=(length(M)+1);
for a=1:length(M)
    set(M(a),'position',[w/c*(a-1) 0.5 w/c 0.5]);
    set(G(a),'position',[w/c*(a-1) 0 w/c 0.5]);
end

% add new ones
N=copyobj(H.magtopo,991);
set(N,'position',[w-w/c 0.5 w/c 0.5],'tag','mags');
% if ~isempty(clim) && isstruct(clim) && isfield(clim,'mags')
%     pos=get(N,'position');
%     C=colorbar; drawnow
%     set(N,'position',pos);
%     colormap(map)
%     set(C,'ytick',(sort([0 clim.mags])-clim.mags(1))/range(clim.mags)*64,'xtick',[]);
%     set(C,'yticklabel',sort([0 clim.mags]));
% end

if meg
    N=copyobj(H.gradtopo,991);
    set(N,'position',[w-w/c 0 w/c 0.5],'tag','grads');
    % if ~isempty(clim) && isstruct(clim) && isfield(clim,'grads')
    %     pos=get(N,'position');
    %     C=colorbar('peer',N); drawnow
    %     set(N,'position',pos);
    %     colormap(hot)
    %     set(C,'ytick',(sort([0 clim.mags])-clim.mags(1))/range(clim.mags)*64,'xtick',[]);
    %     set(C,'yticklabel',sort([0 clim.grads]));
    % end
end

return

%%
function copywave(D)
global H meg magname gradname e
% prevent this if xaxes are inconsistent?
% should colour code by event and indicate sensors by style?
% Legends need some work!
colors=[0         0    1.0000;
    0    0.5000         0;
    1.0000    0         0;
    0    0.7500    0.7500;
    0.7500    0    0.7500;
    0.7500    0.7500    0;
    0.2500    0.2500  0.2500
    0         0         0];
styles={ ...
    {'linestyle','-','linewidth',2,'marker','none'}, ...
    {'linestyle','-','linewidth',1,'marker','none'}, ...
    {'linestyle','-','linewidth',1,'marker','.'}, ...
    {'linestyle',':','linewidth',1,'marker','none'}, ...
    };
% get selected waves
M=findobj(H.selection,'linestyle','-','color',[0 0 0.7]);
G=findobj(H.selection,'linestyle','-','color','k');

figure(992)
Z=findobj(992,'tag','zeroline');
L=findobj(992,'tag','mylegend');
if exist('Z','var') && all(ishandle(Z)); delete(Z); end
if exist('L','var') && all(ishandle(L)); delete(L); end

subplot 211; hold on
U=get(gca,'userData');
N=copyobj(M,gca);
if ~meg; set(gca,'Ydir','reverse'); end %% for compatibility with Vogel
if ~isstruct(U); % first time: set up axes
    U=struct('e',[],'m','');
    ylabel(D.units);
    set(gca,'xtick',get(H.gfp,'xtick'),'xtickLabel',get(H.gfp,'xticklabel'));
    % else
    %     U.e=unique([U.e e]);
    %     U.m=unique([U.m {magname}]);
end
oe=find(ismember(U.e,e));
om=find(ismember(U.m,{strcat(magname{:})}));
if ~isempty(oe); set(N,'color',colors(oe,:));
else
    set(N,'color',colors(length(U.e)+1,:));
    U.e=[U.e e]; U.eh=N;
end
if ~isempty(om); set(N,styles{om}{:});
else
    set(N,styles{length(U.m)+1}{:});
    U.m=[U.m {strcat(magname{:})}]; U.mh=N;
end
set(gca,'userData',U);
% L=legend(U.eh,cellfun(@num2str,num2cell(U.e),'UniformOutput',false));
% set(L,'tag','mylegend');
% L2=legend(U.mh,U.m,'location','northwest');
% set(L2,'tag','mylegend');
axis tight;
Z=plot(repmat(D.events.start,[1 2]),ylim,'k--'); set(Z,'tag','zeroline');
Z=plot(xlim,[0 0],'k--'); set(Z,'tag','zeroline');

subplot 212; hold on
U=get(gca,'userData');
if ~isstruct(U);
    U=struct('e',[],'g','');
    if meg; ylabel('Gradiometer magnitude')
    else ylabel('Contra-ipsi')
    end
    set(gca,'xtick',get(H.gfp,'xtick'),'xtickLabel',get(H.gfp,'xticklabel'));
    %     else
    %         U.e=unique([U.e e]);
    %         U.g=unique([U.g {gradname}]);
end
N=copyobj(G,gca);
oe=find(ismember(U.e,e));
if meg; og=find(ismember(U.g,{strcat(gradname{:})}));
else og=om;
end
if ~isempty(oe); set(N,'color',colors(oe,:));
else
    set(N,'color',colors(length(U.e)+1,:));
    U.e=[U.e e]; U.eh=N;
end
if ~isempty(og); set(N,styles{og}{:});
else
    set(N,styles{length(U.g)+1}{:});
    if meg; U.g=[U.g {strcat(gradname{:})}];
    else U.g=[U.g {strcat(magname{:})}];
    end
    U.gh=N;
end
set(gca,'userData',U);
%     L=legend(U.eh,cellfun(@num2str,num2cell(U.e),'UniformOutput',false));
%     set(L,'tag','mylegend');
%     L2=legend(U.gh,U.g,'location','northwest');
%     set(L2,'tag','mylegend');
axis tight
if ~meg; set(gca,'Ydir','reverse'); end %% for compatibility with Vogel
Z=plot(repmat(D.events.start,[1 2]),ylim,'k--'); set(Z,'tag','zeroline');
Z=plot(xlim,[0 0],'k--'); set(Z,'tag','zeroline')
return

%%
function A = lab2ann(hAx,Hobj,A)
% transfer axis label to parent figure
e=get(Hobj,'extent');
%% Get limits
axun = get(hAx,'Units');
set(hAx,'Units','normalized');  % Need normaized units to do the xform
axpos = get(hAx,'Position');
axlim = axis(hAx);              % Get the axis limits [xlim ylim (zlim)]
axwidth = diff(axlim(1:2));
axheight = diff(axlim(3:4));
%% Transform data from figure space to data space
%x= (x-axlim(1))*axpos(3)/axwidth + axpos(1);
%y= (y-axlim(3))*axpos(4)/axheight + axpos(2);

% Transform and return a position rectangle
pos(1) = (e(1)-axlim(1))/axwidth*axpos(3) + axpos(1);
pos(2) = (e(2)-axlim(3))/axheight*axpos(4) + axpos(2);
pos(3) = e(3)*axpos(3)/axwidth;
pos(4) = e(4)*axpos(4)/axheight;
%% Restore axes units
set(hAx,'Units',axun)
%% create annotation
if pos(1)+pos(3)>1; pos(1)=1-pos(3);
elseif pos(1)<0; pos(1)=0;
elseif pos(1)<0.5; pos(1)=pos(1)-0.004;
end
A(end+1)=annotation('textbox',pos,'string',get(Hobj,'string'),'edge','none');
if get(Hobj,'rotation')==90
    C=get(A(end),'Children');
    set(C(1),'rotation',90)
    set(A(end),'verticalAlignment','middle','horizontalAlignment','center')
end
delete(Hobj)
return