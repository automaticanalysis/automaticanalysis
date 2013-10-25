function RealPar = aas_plot_realignPars(aap,subj,sess,toDisp)
% Based on SPM Realign

% Threshold for excessive movement
QA_TRANSL = 2;
QA_ROT = 8;

EPIs = aas_getimages_bystream(aap,subj,sess,'epi');

[pth, nme, ext] = fileparts(EPIs(1,:));
Params = load(fullfile(pth, ['rp_' nme '.txt']));

% THIS DOES NOT WORK, FLAT PARS... [AVG]
%{
try
    P = spm_vol(EPIs);
    
    if length(P)<2, return; end;
    Params = zeros(numel(P),12);
    for p=1:numel(P),
        Params(p,:) = spm_imatrix(P(subj).mat/P(1).mat);
    end
catch
    [pth, nme, ext] = fileparts(EPIs(1,:));
    Params = load(fullfile(pth, ['rp_' nme '.txt']));
end
%}
RealPar = horzcat(Params(:,1:3),Params(:,4:6)*180/pi);
if toDisp
    if (isempty(spm_figure('FindWin')))
        spm('fmri');
    end;
    fg=spm_figure('FindWin','Graphics');
    if ~isempty(fg)
        % display results
        %-------------------------------------------------------------------
        spm_figure('Clear','Graphics');
        
        % translation over time series
        ax=axes('Position',[0.1 0.65 0.8 0.2],'Parent',fg,'XGrid','on','YGrid','on');
        plot(Params(:,1:3),'Parent',ax); hold on;
        plot(1:size(Params,1),-QA_TRANSL*ones(1,size(Params,1)),'k','Parent',ax)
        plot(1:size(Params,1),QA_TRANSL*ones(1,size(Params,1)),'k','Parent',ax)
        s = ['x translation';'y translation';'z translation'];
        legend(ax, s, 0, 'Location','NorthWest')
        set(get(ax,'Title'),'String','translation','FontSize',16,'FontWeight','Bold');
        set(get(ax,'Xlabel'),'String','image');
        set(get(ax,'Ylabel'),'String','mm');
        YL = get(ax,'YLim');
        ylim([min(-(QA_TRANSL+0.5),YL(1)) max((QA_TRANSL+0.5),YL(2))]);
        
        % rotation over time series
        ax=axes('Position',[0.1 0.35 0.8 0.2],'Parent',fg,'XGrid','on','YGrid','on');
        plot(Params(:,4:6)*180/pi,'Parent',ax); hold on;
        plot(1:size(Params,1),-QA_ROT*ones(1,size(Params,1)),'k','Parent',ax)
        plot(1:size(Params,1),QA_ROT*ones(1,size(Params,1)),'k','Parent',ax)
        s = ['pitch';'roll ';'yaw  '];
        legend(ax, s, 0, 'Location','NorthWest')
        set(get(ax,'Title'),'String','rotation','FontSize',16,'FontWeight','Bold');
        set(get(ax,'Xlabel'),'String','image');
        set(get(ax,'Ylabel'),'String','degrees');
        YL = get(ax,'YLim');
        ylim([min(-(QA_ROT+0.5),YL(1)) max((QA_ROT+0.5),YL(2))]);
        
        % scan-to-scan displacement over time series
        ax=axes('Position',[0.1 0.05 0.8 0.2],'Parent',fg,'XGrid','on','YGrid','on');
        % scale rotation to translation based on the ratio of thresholds (see above)
        plot(diff(sum(abs(horzcat(RealPar(:,1:3),RealPar(:,4:6)/(QA_ROT/QA_TRANSL))),2)),'Parent',ax); hold on;
        plot(1:size(Params,1),-QA_TRANSL*ones(1,size(Params,1)),'k','Parent',ax)
        plot(1:size(Params,1),QA_TRANSL*ones(1,size(Params,1)),'k','Parent',ax)
        set(get(ax,'Title'),'String','Scan-to-scan displacement','FontSize',16,'FontWeight','Bold');
        set(get(ax,'Xlabel'),'String','image');
        set(get(ax,'Ylabel'),'String','a.u. (mm + scaled degrees)');
        YL = get(ax,'YLim');
        ylim([min(-(QA_TRANSL+0.5),YL(1)) max((QA_TRANSL+0.5),YL(2))]);
        
        print('-djpeg','-r75',fullfile(aas_getsubjpath(aap,subj),...
            ['diagnostic_aamod_realign_' aap.acq_details.sessions(sess).name '.jpg']));
    end
end