function fig = meeg_diagnostics_CF(data,diag,savepath)

if any(strcmp(strsplit(data.dimord,'_'),'phase'))
    aas_log([],false,'WARNING: visualising phase-resolved CF data is not yet implemented')
    return
end

multicfg.parameter = 'crsspctrm';
multicfg.showlabels = 'yes';
multicfg.showoutline = 'yes';
multicfg.interactive = 'no';

% save time or create if none
if isfield(data,'time')
    origTime = data.time;
    if isempty(diag.snapshottwoi), diag.snapshottwoi = origTime([1 end]); end
else
    origTime = 0;
    diag.snapshottwoi = [0 0];
end
diag.snapshottwoi = diag.snapshottwoi./1000; % second

% tweak the data into TFR
cfp = data;
cfp.dimord = 'chan_freq_time';
cfp.time = cfp.freqlow;
cfp.freq = cfp.freqhigh;
cfp.crsspctrm = permute(cfp.crsspctrm,[1 3 2 4]);
cfp = removefields(cfp,{'freqlow','freqhigh'});

figName = {''};
if isfield(data,'labelcmb')
    dat = cfp; dat.label = {''}; clear cfp;
    figName = unique(dat.labelcmb(:,1))';
    cfp(1:numel(figName)) = dat;
    for ch = 1:numel(figName)
        ind = strcmp(cfp(ch).labelcmb(:,1),figName{ch});
        cfp(ch).label = union(cfp(ch).labelcmb(ind,2),getedgeelec(cfp(ch).elec),'stable'); % ensure proper arrangement of channels
        cfp(ch).crsspctrm(~ind,:,:,:) = [];
        cfp(ch).crsspctrm(end+1:numel(cfp(ch).label),:,:,:) = NaN;
    end
    cfp = removefields(cfp,'labelcmb');
end

for p = 1:numel(cfp)
    nPlot = size(diag.snapshottwoi,1);
    fig(p) = figure;
    if nPlot > 1
        warning('Figure is too large to fit to screen -> It will be hidden.')
        set(fig(p),'Visible','off');
    end
    set(fig(p),'Position',[0,0,720 nPlot*720]);

    for t = 1:nPlot
        cfg = multicfg;
        cfg.figure = subplot(nPlot,1,t);
        cfg.channel = cfp(p).label;
        tind = (origTime >= diag.snapshottwoi(t,1)) & (origTime <= diag.snapshottwoi(t,2));
        if ~any(tind), continue; end
        tmpcf = cfp(p); tmpcf.crsspctrm = mean(tmpcf.crsspctrm(:,:,:,tind),4);
        ft_multiplotTFR(cfg,tmpcf);
    end
    
    set(fig(p),'Name',figName{p})
    
    if nargin == 3 && ~isempty(savepath)
        figFn = savepath;
        if ~isempty(figName{p}), figFn = [figFn '_' figName{p}]; end
        set(fig(p),'PaperPositionMode','auto');
        print(fig(p),'-noui',[figFn '_multiplot.jpg'],'-djpeg','-r300');
        close(fig(p));
    end
end
end

function ee = getedgeelec(elec)
e1 = min(elec.elecpos);
e2 = max(elec.elecpos);
ee = {...
    elec.label{elec.elecpos(:,2) == e2(2)} ...
    elec.label{elec.elecpos(:,2) == e1(2)} ...
    elec.label{elec.elecpos(:,1) == e2(1)} ...
    elec.label{elec.elecpos(:,1) == e1(1)} ...
    };
end