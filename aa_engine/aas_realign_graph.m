function aas_realign_graph(mv)
% mv can be an array (strvcat or cell) of file names of movement parameters...
% ... or a cell array of matrices containing the parameters...

%% Threshold for excessive movement
QA_TRANSL = 2;
QA_ROT = 8;

%% Input
if ~iscell(mv)
    mv = strvcat2cell(mv);
end

movePars = [];
for s = 1:length(mv)
    if isstr(mv{s})
        % Typically, movePars will be a text file...
        tmp = load(deblank(mv{s}));
        movePars = [movePars; tmp];
    else
        movePars = mv{s};
    end
end

%% Data
% Movements
movePars=movePars-repmat(movePars(1,:),[size(movePars,1) 1]);

movePars(:,4:6)=movePars(:,4:6)*180/pi; % convert to degrees!

mvmean=mean(movePars);
mvmax=max(movePars);
mvstd=std(movePars);

DmovePars = movePars(:,1:3);
RmovePars = movePars(:,4:6);

% scan-to-scan displacement over time series
StS = diff(sum(abs(horzcat(DmovePars,RmovePars/(QA_ROT/QA_TRANSL))),2));

%% Plot
fg= spm_figure;

% translation over time series
ax=axes('Position',[0.1 0.65 0.8 0.2],'Parent',fg,'XGrid','on','YGrid','on');
plot(DmovePars,'Parent',ax); hold on;
plot(1:size(DmovePars,1),-QA_TRANSL*ones(1,size(DmovePars,1)),'k','Parent',ax)
plot(1:size(DmovePars,1),QA_TRANSL*ones(1,size(DmovePars,1)),'k','Parent',ax)
set(get(ax,'Title'),'String','Displacement (mm) [x: blue; y: green; z:red]','FontSize',16,'FontWeight','Bold');
set(get(ax,'Xlabel'),'String','image');
set(get(ax,'Ylabel'),'String','mm');
YL = get(ax,'YLim');
ylim([min(-(QA_TRANSL+0.5),YL(1)) max((QA_TRANSL+0.5),YL(2))]);

% rotation over time series
ax=axes('Position',[0.1 0.35 0.8 0.2],'Parent',fg,'XGrid','on','YGrid','on');
plot(RmovePars*180/pi,'Parent',ax); hold on;
plot(1:size(RmovePars,1),-QA_ROT*ones(1,size(RmovePars,1)),'k','Parent',ax)
plot(1:size(RmovePars,1),QA_ROT*ones(1,size(RmovePars,1)),'k','Parent',ax)
set(get(ax,'Title'),'String','Rotation (deg) [r: blue; p: green; j:red]','FontSize',16,'FontWeight','Bold');
set(get(ax,'Xlabel'),'String','image');
set(get(ax,'Ylabel'),'String','degrees');
YL = get(ax,'YLim');
ylim([min(-(QA_ROT+0.5),YL(1)) max((QA_ROT+0.5),YL(2))]);

% scan-to-scan displacement over time series
ax=axes('Position',[0.1 0.05 0.8 0.2],'Parent',fg,'XGrid','on','YGrid','on');
% scale rotation to translation based on the ratio of thresholds (see above)
plot(StS,'Parent',ax); hold on;
plot(1:size(StS,1),-QA_TRANSL*ones(1,size(StS,1)),'k','Parent',ax)
plot(1:size(StS,1),QA_TRANSL*ones(1,size(StS,1)),'k','Parent',ax)
set(get(ax,'Title'),'String','Scan-to-scan displacement (a.u.)','FontSize',16,'FontWeight','Bold');
set(get(ax,'Xlabel'),'String','image');
set(get(ax,'Ylabel'),'String','a.u. (mm + scaled degrees)');
YL = get(ax,'YLim');
ylim([min(-(QA_TRANSL+0.5),YL(1)) max((QA_TRANSL+0.5),YL(2))]);
    
    
       