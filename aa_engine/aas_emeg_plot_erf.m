function aas_emeg_plot_erf(varargin)
% flip polarity rather than rectify
% !!! for lateralised responses, is it appropriate to directly subtract right
% sensors from left sensors?...
% Assume dipoles will be laterally symmetric;
% then signal at contralateral magnetometers and longitudinal gradiometers
% should have opposite polarity, so should probably flip these;
% contralateral latitudinal gradiometers should have equal polarity already.
% Olaf suggested rectifying the signals prior to subtraction.
% Lateralised analysis perhaps more appropriate in source space.
% In sensor space it would be best if MaxMove has been used to position the
% virtual head symmetrically within the helmet.
% Also, RMS may be best option for averaging over multiple sensors?
%
% Danny Mitchell 2008

addpath /imaging/local/spm/spm5
addpath /imaging/local/spm/spm5/cbu_updates/

if length(varargin)==1; D=varargin{1};
else D=spm_select(1,'mat','Please select a file containing averaged waveforms','',pwd,'^.*m.*mat$');
end

if ~isstruct(D); D=spm_eeg_ldata(char(D)); end

try D.events.repl; if min(D.events.repl)==1; catchnow; end
catch fprintf('\nFound no averaged events. Nothing will be plotted.'); return
end

if D.Nevents>16; fprintf('\nFound >16 events in this file. Nothing will be plotted'); return; end

warning off all; delete(fullfile(D.path,[D.fname(1:end-4) '_sensors.ps'])); warning on all;

% get LHS and RHS which index mags, then longitudinal grads then
% latitudinal grads, in the same order for each hemisphere
[LHS, RHS, trueloc]=getLRpairs(D);

% load data if necessary
if isempty(D.data); D=spm_eeg_ldata(fullfile(D.path,D.fname)); end

% preload file array
D.data=D.data(:,:,:);

% for each event:
for e=1:D.Nevents
    % plot timecourse for:
    % rows: all mags; all grads; peak mag; peak grad; RMS over selected mags; RMS over selected grads
    % columns: both hemisheres; left sensors-right sensors
    warning off all; Fgraph  = spm_figure('GetWin','Graphics'); spm_figure('Clear',Fgraph); warning on all
    set(Fgraph,'renderer','painters','menuBar','none');
    for p=1:12
        leg='';
        subplot(6,2,p)
        data=squeeze(D.data(:,:,e));
        %eog=data(307:308,:);
        if mod(p-1,4)>1 % discard grads and EOG from plots 3,4,7,8,11,12
            data(103:end,:)=0;
            txt='Magnetometer';
        else % discard mags and EOG from other plots
            data(1:102,:)=0;
            data(307:308,:)=0;
            txt='Gradiometer';
        end
        try channames=D.channels.name;
        catch; keyboard % wassup?
        end

        % ldata=abs(ldata); rdata=abs(rdata); % could rectify here rather than trying to account for polarity
        % ...or flip polarity of right-hand magnetometers and longitudinal gradiometers...
        % rdata(1:96,:)=rdata(1:96,:)*-1;
        % data(RHS(1:96,:))=data(RHS(1:96,:))*-1;
        % just do this for right hand graphs

        ldata=data(LHS,:); rdata=data(RHS,:);
        if mod(p,2)==0 % in right hand graphs...
            data=[ldata(1:96,:)+rdata(1:96,:); ldata(97:end,:)-rdata(97:end,:)]; % ... subtract right from left sensors, with flip
            txt=['L-R: ' txt];
            channames=[char(channames(LHS)) repmat('-',length(LHS),1) char(channames(RHS))];
        else channames=char(channames);
        end
        if p>4 % just plot selected sensors for lower 8 graphs
            linewidth=1;
            if p<9 % select sensor with peak range
                range=max(data,[],2)-min(data,[],2);
                [x, ind]=max(range);
                data=data(ind,:);
                txt={[txt ' with peak range: '], channames(ind,:)};
            else % average over all sensors or those closest to coordinates specified in D.ROIs
                ind=any(data')';              
                if isfield(D,'ROIs') %e.g. D.ROIs.space='MNI'; D.ROIs.co=[45, -20, 55];
                    % get coordinates in MNI space (mm)
                    if strcmpi(D.ROIs.space,'MNI'); mni=D.ROIs.co;
                    else mni=tal2mni(D.ROIs.co);
                    end
                    try % convert to subject's MRI space if provided
                        if ~isempty(D.inv{1}.mesh.def);
                            mri=spm_get_orig_coord(mni,D.inv{1}.mesh.def);
                        else mri=mni;
                        end
                        % convert to MEG head coordinates
                        meg=inv(D.inv{1}.datareg.eeg2mri)*[mri'; ones(1,size(mri,1))];
                    catch % Head space origin on MNI152 is approx 0,-16,-48:
                        trans=[1 0 0 0;
                            0 1 0 16;
                            0 0 1 48;
                            0 0 0 1];
                        meg=trans*[mni'; ones(1,size(mni,1))];
                    end
                    % find 12 closest sensors (in head coordinates) to this
                    % coordinate on left
                    [l lord]=sort((trueloc(1,LHS)+abs(meg(1))).^2+(trueloc(2,LHS)-meg(2)).^2+(trueloc(3,LHS)-meg(3)).^2);
                    ind=lord(1:12);
                    txt={txt,sprintf('closest to +/-%g,%g,%g(%s)',abs(D.ROIs.co(1)),D.ROIs.co(2),D.ROIs.co(3),D.ROIs.space{1})};
                else txt={txt};
                end
                data=data(ind,:);
                data(~any(data'),:)=[]; % remove any rows of zeros
                txt{1}=['RMS over selected ' txt{1} 's'];
                data=sqrt(mean(data.^2,1)); % RMS
            end
            if mod(p,2)==0 % include (possibly flipped*) left and right sensors in right-hand graphs
                ldata=ldata(ind,:); rdata=rdata(ind,:);
                ldata(~any(ldata'),:)=[]; % remove any rows of zeros
                rdata(~any(rdata'),:)=[];
                if p>8
                    ldata=sqrt(mean(ldata.^2,1)); % RMS
                    rdata=sqrt(mean(rdata.^2,1)); % RMS
                end
                data=[ldata; data; rdata];
                leg={'L','L-R*','R'};
                ann=annotation('textbox',[0,0,1,0.02],'string','(*Sign inverted for magnetometers and longitudinal gradiometers on right)','edge','none');
            else
                if p==7 % peak magnetometer
                    % plot EOG (is there a better place?) e.g...
                    % subplot('position',[0.6 0.03 0.38 0.06])
                    % subplot('position',[0.57 0.03 0.37 0.06])
                    % data=[data; eog1*range(eog1)/std(eog1)]; % better scaling?
                    % CI would be nice...
                end
            end
        else linewidth=0.5; txt=[txt 's'];
        end

        % modify position of axes
        pos=get(gca,'Position');
        pos(1)=(pos(1)-0.075)*1.1;
        pos(2)=(pos(2)-0.05)*1.05;
        pos(3)=pos(3)*1.2 ;
        set(gca,'Position',pos);

        % when plotting selected sensors (lower 8), superimpose axial plot of sensor position
        if p>4
            plot(trueloc(1,:),trueloc(2,:),'.','Color',[.8 .8 .8]); axis image off; hold on;
            if p==5||p==7||(p==9 && length(ind)>length(LHS))||(p==11 && length(ind)>length(LHS))
                plot(trueloc(1,ind),trueloc(2,ind),'o','color',[0.3 0.3 0.3]);
            else
                leftloc=trueloc(:,LHS(ind)); rightloc=trueloc(:,RHS(ind));
                if mod(p,2)==0           
                    plot(leftloc(1,:),leftloc(2,:),'o','color',[0.4 0.4 1]); %plot(selection(1,:),selection(2,:),'.','color',[0.5 0.5 1]);
                    plot(rightloc(1,:),rightloc(2,:),'o','color',[1 0.4 0.4]); %plot(selection(1,:),selection(2,:),'.','color',[1 0.5 0.5]);
                else
                    plot([leftloc(1,:) rightloc(1,:)],[leftloc(2,:) rightloc(2,:)],'o','color',[0.2 0.2 0.2]); %plot(selection(1,:),selection(2,:),'.','color',[0.5 0.5 1]);
                end                    
            end
            axes('position',get(gca,'Position'));
        end

        % plot data
        data(~any(data'),:)=[]; % remove any rows of zeros
        try 
            if p==1 || p==3
                plot(ldata','linewidth',linewidth,'color','b'); hold on
                plot(rdata','linewidth',linewidth,'color','r');               
            else lines=plot(data','linewidth',linewidth);
                if p==2||p==4; set(lines,'color',[0 0.75 0]); 
                elseif ~mod(p,2)==0; set(lines,'color',[0 0 0]);
                end
            end
        catch % no data ???
        end
        % format and add title
        axis tight; box off; set(gca,'color','none','FontSize',8); title(txt)
        
        % add ticks and time units if final row
        xlims=xlim;
        set(gca,'xtick',xlims(1):(0.2*D.Radc):xlims(2),'xticklabel',[],'fontname','times','tickdir','out');
        if p==12
            if ~isempty(leg);
                legH=legend(leg,'orientation','horizontal','fontweight','bold');
                legend boxoff;
                set(legH,'position',[0.55 0 0.4 0.08]);
            end
        elseif p>10, xlabel('Time (ms)')
            set(gca,'xticklabel',(-D.events.start: 0.2*D.Radc :D.events.stop)*1000/D.Radc);
        end

        % add zero lines
        hold on; plot(ylim,[0 0],'k:');
        plot([D.events.start D.events.start],ylim,'k:');
        drawnow
    end

    %     % could convert to image to save space?
    %     f=getframe(gcf);
    %     subplot('position',[0 0 1 1]);
    %     image(f.cdata); axis off

%% label page
    try labs=D.events.channames{e}; catch; labs=''; end
    try repltxt=sprintf('(%g replications)',D.events.repl(e)); catch; repltxt=''; end
    anntext={sprintf('File: %s',fullfile(D.path,D.fname)), ...
        sprintf('Event/contrast: %g: %s %s.    Date: %s',e,labs,repltxt,datestr(now,0))};
    ann=annotation('textbox',[0.01 0.04  0.91  0.96],'String',anntext,'edge','none');
    
%% print
    [pth nam]=fileparts(D.fname);
    if exist(fullfile(D.path,'figures'),'dir')
        outplot=fullfile(D.path,'figures',[nam '_sensors.ps']);
    else outplot=fullfile(D.path,[nam '_sensors.ps']);
    end

    spm_print(outplot); % does something to avoid segmentation fault?
    print('-dpng', '-zbuffer','-r82', regexprep(outplot,'.ps',sprintf('_%02g.png',e))); % low res, but filts in web browser
    delete(ann)
end

return
