function aas_emeg_plot_topographies(varargin)
% plot_topographies()
% plot_topographies(D,movie,T,sensors,figstyle,layout,fulloutputfile,vector)
%
% Plot 2D projection of sensor data at selected time points, averaged
% across time-window, or as a movie. Various layout options. Plot
% either magnetometer signal, or power for gradiometer pairs.
%
% D can be structure or filename or cell array of these.
% movie  - 'off' or frame length in ms
% T - specific time (ms), or
%   - time window (ms) to average, or
%   - vector of snapshots
% sensors - 'Mags' , 'Grads' or 'EEG'
% figstyle - 'A4' or 'Screen'
% layout - 'subs x time'
%        - 'time x subs'
%        - 'subs x evts'
%        - 'evts x subs'
%        - 'evts x time'
%        - 'time x evts'
%        - 'subs x subs'
%        - 'time x time'
%        - 'evts x evts'
% fulloutputfile - full path of output file (ending will be
% changed to <sensors>_<event>.png/.ps for figure,
% or .avi for movie)
% vector - 0|1 flag to plot as image (bitmap) or surface (vector) respectively
%
% Danny Mitchell 20/01/08

% Could confirm files are correct
% Improve marking of sensor positions
% Would be nice to add contour plot of stats for group analysis
% Using cubic interpolation which is smoother and faster than default 
% linear interpolation, but seems to have lower contrast, probably due to smoothing.
% Use all sensors for mags or grads, and set bad eeg channels to zero.

load /imaging/dm01/MEG/aaMEG/redblue2_128
addpath /imaging/local/spm/spm5
addpath /imaging/local/spm/spm5/cbu_updates/

if length(varargin)>0; data=varargin{1};
else data=spm_select(Inf,'mat','Please select file(s) containing averaged waveforms','',pwd,'^.*m.*mat$');
end
if ischar(data(1,:)); D=spm_eeg_ldata(deblank(char(data(1,:))));
else D=data; data=fullfile(D.path,D.fname);
end

if length(varargin)>1; framelength=varargin{2};
else framelength=spm_input('Create movie?',1,'bn1','None','off',4);
end

if ischar(framelength);
    if length(varargin)>2; T=varargin{3};
    elseif isempty(varargin); T=spm_input('Times(s) or timewindow (ms):',1,'e',0);
    end
else
    if length(varargin)>2; T=varargin{3};
    elseif isempty(varargin); T=spm_input('Times or timewindow (ms) from which to set colour bar:',1,'e',0);
    end
end

Tms=T;
T = round(Tms * D.Radc/1000) + D.events.start+1; % samples

if length(T)==2; T={T(1):T(2)}; timelabel=sprintf('Mean of %g to %g ms',Tms(1),Tms(2));
else timelabel=sprintf('%g ms',Tms(1));
end

if length(varargin)>3; sensors=varargin{4};
elseif isempty(varargin); sensors=spm_input('Sensor type:','+1','b','Mags|Grads|EEG');
end

if length(varargin)>4; figstyle=varargin{5};
else figstyle=char(spm_input('Figure style:','+1','b','A4|Screen',{'A4','Screen'},1));
end

layouts={'subs x time' ... % 1
    'time x subs' ... % 2
    'subs x evts' ... % 3
    'evts x subs' ... % 4
    'evts x time' ... % 5
    'time x evts' ... % 6
    'subs x subs' ... % 7
    'time x time' ... % 8
    'evts x evts'}; % 9
if isnumeric(framelength);
    layouts([1 2 5 6 8])=[];
    timelabel=['Colour scale estimated from ' timelabel];
end

if length(varargin)>5; layout=varargin{6};
else layout=char(spm_input('Layout:','+1','m',layouts,layouts,length(layouts)));
end

[pa na]=fileparts(D.fname);
if length(varargin)>6; fulloutfile=varargin{7};
else fulloutfile=spm_input('Full output path:','+1','s',[na '_topo_' sensors '.png']);
end

if length(varargin)>7; vector=varargin{8};
else vector=0;
end

%%%%%%%% following adapted from spm_eeg_scalp2d_ext
CTF = load(fullfile(spm('dir'), 'EEGtemplates', D.channels.ctf));
if strcmp(sensors,'EEG'); 
    Gfx.channels=setdiff(1:length(D.channels.eeg),D.channels.Bad);
else % just use 1st 102 MEG channels
    Gfx.channels=1:102;
end

CTF.Cpos = CTF.Cpos(:, D.channels.order(Gfx.channels));

handles.x = min(CTF.Cpos(1,:)):0.005:max(CTF.Cpos(1,:)); % need to be double for griddata
handles.y = min(CTF.Cpos(2,:)):0.005:max(CTF.Cpos(2,:));
[handles.x1, handles.y1] = meshgrid(handles.x, handles.y);
handles.xp = CTF.Cpos(1,:)';
handles.yp = CTF.Cpos(2,:)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=1; t=1; e=1;
switch layout
    case 'subs x subs'
        nrows = round(size(data,1).^0.5);
        ncols = ceil(size(data,1).^0.5);
        npags = length(T)*D.Nevents;
        t=0;e=1;
    case 'time x time'
        nrows = round(length(T).^0.5);
        ncols = ceil(length(T).^0.5);
        npags = size(data,1)*D.Nevents;
        e=0;s=1;
    case 'evts x evts'
        nrows = round(D.Nevents.^0.5);
        ncols = ceil(D.Nevents.^0.5);
        npags = length(T)*size(data,1);
        s=1;t=0;
    otherwise
        switch layout(1:4)
            case 'subs'; nrows=size(data,1);
            case 'time'; nrows=length(T);
            case 'evts'; nrows=D.Nevents;
        end
        switch layout(8:11)
            case 'subs'; ncols=size(data,1);
            case 'time'; ncols=length(T);
            case 'evts'; ncols=D.Nevents;
        end
        npags=1;
        if isempty(strfind(layout,'subs')); npags=npags*size(data,1); end
        if isempty(strfind(layout,'time')); npags=npags*length(T); end
        if isempty(strfind(layout,'evts')); npags=npags*D.Nevents; end
end

if strcmp(layout(1:4),layout(8:11)) && strcmp(figstyle,'A4')
    tmp=nrows;
    nrows=ncols;
    ncols=tmp;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nPlotting topographies...'); tic;
for p=1:npags
    pos=0;
    if strcmp(layout,'subs x subs');
        t=t+1;s=0;
        if t>length(T); t=1; e=e+1; end
        try ename=D.events.names{e}; catch ename=num2str(e); end
        header={['Time = ' timelabel], ...
            ['Event = ' ename]};
    elseif strcmp(layout,'time x time');
        e=e+1;t=0;
        if e>D.Nevents; e=1; s=s+1; end
        header={['File = ' data(s,:)], ...
            ['Event = ' ename]};
    elseif strcmp(layout,'evts x evts');
        t=t+1;e=0;
        if t>length(T); t=1; s=s+1; end
        header={['File = ' data(s,:)], ...
            ['Time = ' timelabel]};
    elseif isempty(strfind(layout,'subs')); header={['File = ' data(s,:)]};
    elseif isempty(strfind(layout,'time')); header={['Time = ' timelabel]};
    elseif isempty(strfind(layout,'evts')); try ename=D.events.names{e}; catch ename=num2str(e); end; header={['Event = ' ename]};
    end

    switch figstyle
        case 'A4'
            warning off all;
            Fh=spm_figure('GetWin','Graphics'); spm_figure('Clear',Fh);
            warning on all
            pp=get(Fh,'paperposition');
            set(Fh,'paperposition',[pp(1)/2 pp(2)/2 pp(3)+pp(1) pp(4)+pp(2)],'menuBar','none');
        otherwise
            Fh=figure(10); clf;
            set(0,'Units','pixels');
            screen=get(0,'ScreenSize');
            set(10,'Position',screen,'resize','off');
            clear screen;
    end

    % label page
    header=[header, sprintf('Printed on %s',datestr(now,0))];
    ann=annotation('textbox',[0.01 0.04  0.91  0.96],'String',header,'edge','none','interpreter','none');
    clear header
    
    if strcmp(sensors,'Grads');
        colormap(hot)
    else
    try colormap(map);
    catch colormap jet;
    end
    end
    cmax=-Inf; cmin=Inf; allmax=[];

    for y=1:nrows
        for x=1:ncols
            if strcmp(layout,'subs x subs'); s=s+1;
            elseif strcmp(layout,'time x time'); t=t+1;
            elseif strcmp(layout,'evts x evts'); e=e+1;
            else
                if strcmp(layout(1:4),'subs'); s=y;
                elseif strcmp(layout(8:11),'subs'); s=x;
                else s=p;
                end

                if strcmp(layout(1:4),'time'); t=y;
                elseif strcmp(layout(8:11),'time'); t=x;
                else t=p;
                end

                if strcmp(layout(1:4),'evts'); e=y;
                elseif strcmp(layout(8:11),'evts'); e=x;
                else e=p; try ename=D.events.names{e}; catch ename=num2str(e); end;
                end
            end

            if e>D.Nevents; break; end
            
            pos=pos+1;
            if strcmp(layout(1:4),layout(8:11)) % leave space for titles
                sp{pos}=subplot(nrows,ncols,pos,'xtick',[],'ytick',[]);
            else % tighter format 
                sp{pos}=subplot('position',[0.05+(x-1)*(0.95/ncols) 0.94-y*(0.94/(nrows+1)) 0.95/ncols 0.94/(nrows+1)]);
            end
            
            try
                if ischar(data(s,:)); D=spm_eeg_ldata(deblank(char(data(s,:))));
                else D=data;
                end
            catch
                delete(gca)
                pos=pos-1;
                sp=sp(1:end-1);
                break
            end

            D.data=D.data(:,:,:); % prefetch data
            try % maybe save a little memory
                D=rmfield(D,'filter');
                D=rmfield(D,'ica');
                D=rmfield(D,'thresholds');
                D=rmfield(D,'ROIs');
            catch
            end
            
            if strcmp(sensors,'Grads');
                % (channel sorting from JT):
                name=D.channels.name;
                ch=[];
                for i=1:102;
                    ch=[ch;find(strncmpi(name,name{i}(1:6),6))'];
                end
                 % get pythagorean sum for each gradiometer pair and copy
                 % to 1st 102 channels
                D.data(1:102,:,:)=sqrt(D.data(ch(:,2),:,:).^2+D.data(ch(:,3),:,:).^2);
                clear name ch
            end

            if iscell(T)
                d = squeeze(mean(D.data(Gfx.channels, T{1}, e), 2));
            else
                d = squeeze(D.data(Gfx.channels, T(t), e));
                timelabel=sprintf('%g ms',Tms(t));
            end

            if isnumeric(framelength) && t==1
                alldata{pos} = squeeze(D.data(Gfx.channels,:,e));
            elseif isnumeric(framelength)
                alldata{pos} = NaN;
            end

            if isnan(sum(d(:))); xlabel('No data'); continue; end

            z = griddata(handles.xp, handles.yp, d, handles.x1, handles.y1,'cubic');
            cmax=max([max(z(:)); cmax]);
            cmin=min([min(z(:)); cmin]);
            allmax(end+1)=max(z(:));

            if vector
                surf{pos}=surface(handles.x, handles.y, z,'facelighting','phong');
                axis equal tight off; shading('interp'); hold on
                if strcmp(sensors,'EEG'); % mark electrode positions
                    set(gca,'ycolor',get(Fh,'color'),'xcolor',get(Fh,'color'),'xtick',[],'ytick',[])
                    plot3(handles.xp, handles.yp, d, '.','color',[0.5 0.5 0.5]);
                end
            else
                imagesc(z)
                axis xy equal tight off
            end
            set(sp{pos},'xtick',[],'ytick',[]);

            if strcmp(layout(1:4),layout(8:11))
                if strcmp(layout(1:4),'subs'); title(regexprep(D.path,'^.*/(.*/[^/]*)','$1'),'interpreter','none'); end
                if strcmp(layout(1:4),'time'); title(T(t)); end
                if strcmp(layout(1:4),'evts');
                    try title(D.events.names{e}); catch title(['Event ' num2str(e)]); end
                end
                %set(get(sp{1},'title'),'position',[0.5 0.95 0])
            end

            if x==1 && ~strcmp(layout(1:4),layout(8:11))
                l=ylabel('');
                if strcmp(layout(1:4),'subs'); text('position',get(l,'position'),'string',regexprep(D.path,'^.*/(.*/[^/]*)','$1'),'horizontalalignment','center','rotation',90,'interpreter','none'); end
                if strcmp(layout(1:4),'time'); text('position',get(l,'position'),'string',timelabel,'horizontalalignment','center','rotation',90); end
                if strcmp(layout(1:4),'evts');
                    try text('position',get(l,'position'),'string',D.events.names{e},'horizontalalignment','center','rotation',90,'interpreter','none','fontsize',8);
                    catch text('position',get(l,'position'),'string',['Event ' num2str(e)],'horizontalalignment','center','rotation',90,'interpreter','none','fontsize',8);
                    end
                end
            end

            if y==1 && ~strcmp(layout(1:4),layout(8:11))
                if strcmp(layout(8:11),'subs'); title(regexprep(D.path,'^.*/(.*)/[^/]*','$1'),'interpreter','none'); end
                if strcmp(layout(8:11),'time'); title(timelabel); end
                if strcmp(layout(8:11),'evts');
                    try title(D.events.names{e}); catch title(['Event ' num2str(e)]); end
                end
            end

            drawnow
        end % next column
    end % next row

    %     set(gca,'ycolor','k')
    %     set(gca,'xcolor','k')
    %     set(gca,'ytick',ylim); set(gca,'yticklabel',{'A','P'});
    %     set(gca,'xtick',xlim); set(gca,'xticklabel',{'L','R'});
    %     set(gca,'color','none')
    %     axis on

    % give all axes same colour scale
    if strcmp(sensors,'Grads');
        %myrange=[0 mean(allmax)];
        %myrange=[0 prctile(allmax,200/3)];
        myrange=[0 cmax];
    else % (symmetric about zero)
        m=min(abs([cmin cmax])); myrange=[-m m]; clear m
    end
    if isnumeric(framelength); myrange=myrange*1.5; end % expand caxis for movie since estimated range may be too small?
    for csp=1:length(sp); caxis(sp{csp},myrange); end
    
    cpb=get(sp{ncols},'position');
    gap=max([0.05, 1-cpb(1)-cpb(3)]);

    Hc=colorbar('location','south','position',[0.5+gap/2 cpb(2) 0.5-1.5*gap 0.01],'xaxislocation','bottom');
    if strcmp(sensors,'Grads'); xlabel(Hc, ['(' D.units '/m) ???'], 'FontSize', 12);
    else strcmp(sensors,'Mags'); xlabel(Hc, D.units, 'FontSize', 12);
    end
    op=get(Hc,'outerPosition');set(Hc,'outerPosition',[op(1) 0.03 op(3:4)]);
    clear cpb op
    
    % try colormap(map); catch; end; % for greyscalable
    
    if ~isnumeric(framelength)
        fprintf('printing...')
        % print
        if isnumeric(fulloutfile)
            fulloutfile=fullfile(D.path,sprintf('_%s_topo_%s_%02g.png',na, sensors, e));
        else
            % add sensor type and event and convert extension to png
            part1=regexprep(fulloutfile,{'(_\d+)*\..*$','_?EEG_?','_?Mags_?','_?Grads_?'},{'','','',''});
            if isempty(strfind(part1,'/')); part1=fullfile(D.path,part1); end
            if isempty(strfind(layout,'evts')); fulloutfile=[part1 sprintf('_%s_%02g.png',sensors, e)];
            else fulloutfile=[part1 sprintf('_%s.png',sensors)];
            end
        end
        if strcmp(figstyle,'A4');
            psfile=strrep(fulloutfile,'.png','.ps');
            warning off all; try delete(psfile); catch end; warning on all;
            spm_print_zbuffer(psfile); % does something to avoid segmentation fault?
        end
        print('-dpng','-r82', fulloutfile); % low res, but filts in web browser
    elseif exist('Hc','var')
        % create movie
        if isnumeric(fulloutfile)
            fulloutfile=fullfile(D.path,[na '_topo_' sensors '-movie_' num2str(e) '.avi']);
        else
            % convert extension to png
            fulloutfile=regexprep(fulloutfile,'\..*$','.avi');
        end

        frame=0;
        % setup time marker
        % Htb=subplot('position',[gap cpb(2)-0.02 0.5-1.5*gap 0.01]);
        Htb=subplot('position',[gap 0.03 0.5-1.5*gap 0.01]);
        cbpos=get(Hc,'position'); set(Htb,'Position',[gap cbpos(2:4)]);
        set(Hc,'HandleVisibility','off') % stop it being deleted/overwritten
        clear gap cbpos
        
        tm=plot(-D.events.start,0,'r^','markerfacecolor','r','markersize',10);
        step=ceil((D.events.stop+D.events.start)/5*1000/D.Radc/100)*100*D.Radc/1000;
        ticks=union(handles.ms(1:step:end),0);
        xlim([-D.events.start D.events.stop]);
        set(gca,'xtick',union(0,-D.events.start:step:D.events.stop));
        set(gca,'xticklabel',ticks);
        set(gca,'ytick',[]);
        clear tics step
        
        titles=cell(1,length(sp));
        for csp=1:length(sp)
            titles{csp}=get(get(sp{csp},'title'),'String');
        end
        clear x y tmp
        for sample=1:framelength*D.Radc/1000:D.Nsamples
            frame=frame+1;
            for csp=1:length(sp) % loop through subplot handles
                if isnan(alldata{csp}(1))
                    cla(sp{csp});
                else
                    %try delete(surf{csp}); catch end
                    subplot(sp{csp});
                    z = griddata(handles.xp, handles.yp, alldata{csp}(:,round(sample)), handles.x1, handles.y1,'cubic');
                    if vector
                        set(surf{csp},'Cdata', z);
                    else
                        imagesc(z,myrange)
                        clear z
                        axis xy equal tight off
                        title(titles{csp});
                    end
                end
            end
            set(tm,'xdata',sample-D.events.start); % update time marker
            xlabel(Htb, sprintf('Time = %g ms',(sample-D.events.start)*1000/D.Radc), 'FontSize', 12);
            M(frame)=getframe(Fh); % I think this includes "drawnow"?
        end
        fprintf('\nSaving movie to %s...',fulloutfile)
        warning off all; try delete(fulloutfile); catch end;
        movie2avi(M,fulloutfile);
        warning on all
        clear M

        % compress movie using mplayer??
        compfile=strrep(fulloutfile,'.avi','_C.avi');
        fprintf('\nAttempting to compress AVI using mencoder...')
        cmd=sprintf('/imaging/dm01/MoreTools/linux/MPlayer-1.0rc2/mencoder %s -o %s -ovc lavc -ffourcc DX50',fulloutfile,compfile);
        try unix(cmd);
        catch fprintf('Failed!')
        end
        
        % this used lots of memory and if we try again it will probably 
        % crash catastrophically, so force memory error now...
        error('MATLAB:nomem', ['Forcing this error after creating a video\n' ...
                               'otherwise matlab might crash catastrophically\n' ...
                               'if attempting to create another video without\n' ...
                               'restarting matlab to free up memory.'])
    end

    try delete(ann); catch end;
end % next page

%delete(Fh)
%toc
return