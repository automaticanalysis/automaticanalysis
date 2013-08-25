function aas_emeg_plot_sensor_spms(aap,Overwrite)

if ~exist('aap','var');
    clear all
    [aap ok]=spm_select(1,'any','Please select aap_parameters file or script to generate aap structure','',pwd,'^aa.*\.m.*$');
    if ~ok; return; end
    try load(aap);
    catch error('\nFailed to load aap structure.\n')
    end
end
addpath /imaging/dm01/MEG/aaMEG/

[junk dirs]=spm_select('FPList', ...
    fullfile(aap.acq_details.root,'GroupAnalysis_Sensors'),'');
dirs=cellstr(dirs);

outdir=fullfile(aap.acq_details.root,'GroupAnalysis_Sensors','figures');
if ~exist(outdir,'dir'); mkdir(outdir); end
outdir=fullfile(outdir,'spms');
if ~exist(outdir,'dir'); mkdir(outdir); end

if ~exist('Overwrite','var'), Overwrite=0; end

summaryfile=fullfile(aap.acq_details.root,'GroupAnalysis_Sensors', ...
    'figures','SummaryOfSPMs.ps');
if exist(summaryfile,'file'); delete(summaryfile); end

groups=unique(regexprep(dirs,'-....?.?$','')); % remove sensor suffix

for g=1:length(groups) % ignore . and ..

    warning off all; F=spm_figure('GetWin','Graphics'); spm_figure('Clear',F); warning on all

    if isempty(findstr(groups{g},'2ndLev')); continue; end; % e.g. figures dir

    if ~isempty(findstr(groups{g},'_mt'))
        % assume for now this is tf.
        chantypes=regexprep(dirs(strmatch(groups{g},dirs)),'.*average-','');
        tf=true;
        timewin=regexprep(groups{g},'.*(mt\d).*','$1');
    else
        chantypes=regexprep(dirs(strmatch(groups{3},dirs)),'.*ms-','');
        chantypes(strcmpi('Heog',chantypes))=[];
        chantypes(strcmpi('Veog',chantypes))=[];
        tf=false;
        timewin=regexp(groups{g},'_([^_]*$)','tokens');
        timewin=timewin{1}{1};
    end

    prefixes=regexp(groups{g},'Lev_([^/]*N?ST?t?\d?)_','tokens');
    prefixes=prefixes{1}{1};

    [junk effects]=spm_select('List',regexprep(sprintf('%s-%s',groups{g},chantypes{1}),'-$',''),'');
    effects=cellstr(effects);

    for e=3:length(effects) % ignore . and ..

        %     elseif isempty(findstr(groups{dg},'-eeg'))
        %         load /imaging/local/spm/spm5/EEGtemplates/FIF306_setup.mat
        %         topchan=102;
        %         chantypes={regexprep(groups{d},'.*ms-','')}; % assume filename format
        %         % following chantypes will actually be in their own folders
        %         % chantypes={'mags','grds','mags-RMS','grds-RMS','Heog','Veog'};
        %     else
        %         load /imaging/local/spm/spm5/EEGtemplates/draft_eeg_montage_djm_mean.mat
        %         topchan=size(Cpos,2);
        %         chantypes={'eeg'};
        %     end
        outfile=fullfile(outdir,sprintf('%s_%s.png',timewin,effects{e}));
        if ~exist(outfile,'file') || Overwrite
            fprintf('\nWorking on %s',outfile)
            if tf; clim=[Inf -Inf]; end

            for ct=1:length(chantypes)

                sp{ct}=subplot(ceil(sqrt(length(chantypes))),floor(sqrt(length(chantypes))),ct);
                chantype=chantypes{ct};
                indir=regexprep(sprintf('%s-%s',groups{g},chantype),'-$','');

                if exist(fullfile(indir,effects{e},'beta_0002.img'),'file')
                    Ktype='mean';
                    if exist(fullfile(indir,effects{e},[Ktype 'K.img']),'file') continue; end
                    % factor has more than 2 levels; choose a contrast to display?
                    fprintf('\nComputing K-weighted contrast')
                    txt=sprintf('Superimposed on %s K-weighted contrast',Ktype);
                    % '' for individual K weighted contrasts;
                    % 'mean' to use the mean K function for all subjects

                    try
                        clear K_VSTM* K_ESTA* Y Cnames Cpos
                        % gosh, here comes some ugly code...
                        if ~isempty(regexp(timewin,'saverage','ONCE'))
                            K_VSTM_L=textscan(ls(fullfile(aap.acq_details.root,'*', ...
                                sprintf('FullFact_1stLev_%s*_SampleK_%s*%s*',prefixes,timewin,regexprep(indir,'[^-]*(-.*)','$1')), ...
                                sprintf('%sK_VSTM_L.img',Ktype))),'%s');
                            V=spm_vol(char(K_VSTM_L{1}));
                            K_VSTM_L=mean(single(spm_read_vols(V)),4); clear V

                            K_VSTM_R=textscan(ls(fullfile(aap.acq_details.root,'*', ...
                                sprintf('FullFact_1stLev_%s*_SampleK_%s*%s*',prefixes,timewin,regexprep(indir,'[^-]*(-.*)','$1')), ...
                                sprintf('%sK_VSTM_R.img',Ktype))),'%s');
                            V=spm_vol(char(K_VSTM_R{1}));
                            K_VSTM_R=mean(single(spm_read_vols(V)),4);clear V

                            K_ESTA_L=textscan(ls(fullfile(aap.acq_details.root,'*', ...
                                sprintf('FullFact_1stLev_%s*_SampleK_%s*%s*',prefixes,timewin,regexprep(indir,'[^-]*(-.*)','$1')), ...
                                sprintf('%sK_ESTA_L.img',Ktype))),'%s');
                            V=spm_vol(char(K_ESTA_L{1}));
                            K_ESTA_L=mean(single(spm_read_vols(V)),4);clear V

                            K_ESTA_R=textscan(ls(fullfile(aap.acq_details.root,'*', ...
                                sprintf('FullFact_1stLev_%s*_SampleK_%s*%s*',prefixes,timewin,regexprep(indir,'[^-]*(-.*)','$1')), ...
                                sprintf('%sK_ESTA_R.img',Ktype))),'%s');
                            V=spm_vol(char(K_ESTA_R{1}));
                            K_ESTA_R=mean(single(spm_read_vols(V)),4);clear V
                        elseif ~strcmp(timewin(1),'m')
                            K_VSTM_L=textscan(ls(fullfile(aap.acq_details.root,'*', ...
                                sprintf('FullFact_1stLev_%s*_SampleK_avewin_%s-%s',prefixes,timewin,chantype), ...
                                sprintf('%sK_VSTM_L.img',Ktype))),'%s');
                            K_VSTM_L=spm_read_vols(spm_vol(char(K_VSTM_L{1})));

                            K_VSTM_R=textscan(ls(fullfile(aap.acq_details.root,'*', ...
                                sprintf('FullFact_1stLev_%s*_SampleK_avewin_%s-%s',prefixes,timewin,chantype), ...
                                sprintf('%sK_VSTM_R.img',Ktype))),'%s');
                            K_VSTM_R=spm_read_vols(spm_vol(char(K_VSTM_R{1})));

                            K_ESTA_L=textscan(ls(fullfile(aap.acq_details.root,'*', ...
                                sprintf('FullFact_1stLev_%s*_SampleK_avewin_%s-%s',prefixes,timewin,chantype), ...
                                sprintf('%sK_ESTA_L.img',Ktype))),'%s');
                            K_ESTA_L=spm_read_vols(spm_vol(char(K_ESTA_L{1})));

                            K_ESTA_R=textscan(ls(fullfile(aap.acq_details.root,'*', ...
                                sprintf('FullFact_1stLev_%s*_SampleK_avewin_%s-%s',prefixes,timewin,chantype), ...
                                sprintf('%sK_ESTA_R.img',Ktype))),'%s');
                            K_ESTA_R=spm_read_vols(spm_vol(char(K_ESTA_R{1})));
                        else
                            K_VSTM_L=textscan(ls(fullfile(aap.acq_details.root,'*', ...
                                sprintf('FullFact_1stLev_%s*_SampleK*%s',prefixes,regexprep(indir,'[^-]*(-.*)','$1')), ...
                                sprintf('%sK_VSTM_L.img',Ktype))),'%s');
                            K_VSTM_L=spm_read_vols(spm_vol(char(K_VSTM_L{1})));

                            K_VSTM_R=textscan(ls(fullfile(aap.acq_details.root,'*', ...
                                sprintf('FullFact_1stLev_%s*_SampleK*%s',prefixes,regexprep(indir,'[^-]*(-.*)','$1')), ...
                                sprintf('%sK_VSTM_R.img',Ktype))),'%s');
                            K_VSTM_R=spm_read_vols(spm_vol(char(K_VSTM_R{1})));

                            K_ESTA_L=textscan(ls(fullfile(aap.acq_details.root,'*', ...
                                sprintf('FullFact_1stLev_%s*_SampleK*%s',prefixes,regexprep(indir,'[^-]*(-.*)','$1')), ...
                                sprintf('%sK_ESTA_L.img',Ktype))),'%s');
                            K_ESTA_L=spm_read_vols(spm_vol(char(K_ESTA_L{1})));

                            K_ESTA_R=textscan(ls(fullfile(aap.acq_details.root,'*', ...
                                sprintf('FullFact_1stLev_%s*_SampleK*%s',prefixes,regexprep(indir,'[^-]*(-.*)','$1')), ...
                                sprintf('%sK_ESTA_R.img',Ktype))),'%s');
                            K_ESTA_R=spm_read_vols(spm_vol(char(K_ESTA_R{1})));
                        end

                        if ~isempty(regexp(effects{e},'SepTask_DifCue.*_1','ONCE'))
                            Y=(K_VSTM_L - K_VSTM_R)./2;
                        elseif ~isempty(regexp(effects{e},'SepTask_DifCue.*_2','ONCE'))
                            Y=(K_ESTA_L - K_ESTA_R)./2;
                        elseif ~isempty(regexp(effects{e},'ComTask_DifCue.*','ONCE'))
                            Y=(K_VSTM_L + K_ESTA_L)./2 - (K_VSTM_R + K_ESTA_R)./2;
                        elseif ~isempty(regexp(effects{e},'DifTask_SepCue.*_1','ONCE'))
                            Y=(K_VSTM_L - K_ESTA_L)./2;
                        elseif ~isempty(regexp(effects{e},'DifTask_SepCue.*_2','ONCE'))
                            Y=(K_VSTM_R - K_ESTA_R)./2;
                        elseif ~isempty(regexp(effects{e},'DifTask_ComCue.*','ONCE'))
                            Y=(K_VSTM_L + K_VSTM_R)./2 - (K_ESTA_L + K_ESTA_R)./2;
                        elseif ~isempty(regexp(effects{e},'DifTask_DifCue.*','ONCE'))
                            Y=(K_VSTM_L - K_VSTM_R - K_ESTA_L + K_ESTA_R)./4;
                        elseif ~isempty(regexp(effects{e},'ComTask_ComCue.*','ONCE'))
                            Y=(K_VSTM_L + K_VSTM_R + K_ESTA_L + K_ESTA_R)./4;
                        elseif ~isempty(regexp(effects{e},'SepTask_ComCue.*_1','ONCE'))
                            Y=(K_VSTM_L + K_VSTM_R)./2;
                            %                             %%% debug:
                            %                             nr=ceil(sqrt(size(Y,4)));
                            %                             for p=1:size(Y,4)
                            %                                 subplot(nr,nr,p)
                            %                                 imagesc(squeeze(Y(:,:,1,p)));
                            %                                 lim=max(abs(caxis));
                            %                                 caxis([-lim lim]);
                            %                             end
                            %                             subplot(nr,nr,p+1)
                            %                             imagesc(mean(Y,4));
                            %                             lim=max(abs(caxis));
                            %                             caxis([-lim lim]);
                            %                             xlabel({effects{e},chantype})
                            %                             spm_print(summaryfile);
                            %                             %%%
                        elseif ~isempty(regexp(effects{e},'SepTask_ComCue.*_2','ONCE'))
                            Y=(K_ESTA_L + K_ESTA_R)./2;
                            %                             %%% debug:
                            %                             nr=ceil(sqrt(size(Y,4)));
                            %                             for p=1:size(Y,4)
                            %                                 subplot(nr,nr,p)
                            %                                 imagesc(squeeze(Y(:,:,1,p)));
                            %                                 lim=max(abs(caxis));
                            %                                 caxis([-lim lim]);
                            %                             end
                            %                             subplot(nr,nr,p+1)
                            %                             imagesc(mean(Y,4));
                            %                             lim=max(abs(caxis));
                            %                             caxis([-lim lim]);
                            %                             xlabel({effects{e},chantype})
                            %                             spm_print(summaryfile);
                            %                             %%%
                        elseif ~isempty(regexp(effects{e},'ComTask_SepCue.*_1','ONCE'))
                            Y=(K_VSTM_L + K_ESTA_L)./2;
                        elseif ~isempty(regexp(effects{e},'ComTask_SepCue.*_2','ONCE'))
                            Y=(K_VSTM_R + K_ESTA_R)./2;
                        elseif ~isempty(regexp(effects{e},'SepTask_SepCue.*_1','ONCE'))
                            Y=K_VSTM_L;
                        elseif ~isempty(regexp(effects{e},'SepTask_SepCue.*_2','ONCE'))
                            Y=K_VSTM_R;
                        elseif ~isempty(regexp(effects{e},'SepTask_SepCue.*_3','ONCE'))
                            Y=K_ESTA_L;
                        elseif ~isempty(regexp(effects{e},'SepTask_SepCue.*_4','ONCE'))
                            Y=K_ESTA_R;
                        else
                            debugnow
                            Y=spm_read_vols(spm_vol(fullfile(indir,effects{e},'spmF_0001.img')));
                            txt='Superimposed on F statistic';
                        end

                    catch
                        debugnow
                        Y=spm_read_vols(spm_vol(fullfile(indir,effects{e},'spmF_0001.img')));
                        txt='Superimposed on F statistic';
                    end

                    if ndims(Y)==4
                        Y=mean(Y,4);
                    end
                    
                    if ndims(Y)==2
                        
                        %%% also print the individual topographies
                        tempoutfile=strrep(outfile,'.png',sprintf('_allsubs_%s.ps',chantype)); %%%
                        WS   = spm('WinScale');
                        Rect = spm('WinSize','Graphics','raw').*WS;
                        tempH= figure('Position',Rect,'Resize','off','Color','w','PaperType','A4',...
                            'PaperUnits','normalized','PaperPosition',[.0726 .0644 .854 .870],'Toolbar','none');
                        trows=ceil(sqrt(size(Y,4))); %%%
                        tcols=floor(sqrt(size(Y,4)));
                        if (trows*tcols)<size(Y,4); trows=trows+1; end;
                        tp=0;
                        tclim=[min(Y(:)) max(Y(:))];
                        a=textscan(ls(fullfile(aap.acq_details.root,'*', ...
                            sprintf('FullFact_1stLev_%s*_SampleK*%s*',prefixes,regexprep(indir,'[^-]*(-.*)','$1')), ...
                            sprintf('%sK_VSTM_L.img',Ktype))),'%s');
                        nams=regexprep(a{1},'.*/([^/]*)/FullFact.*','$1');
                        for tsub=1:size(Y,4)
                            tp=tp+1; subplot(trows,tcols,tp)
                            imagesc(squeeze(Y(:,:,:,tsub))); caxis(tclim);
                            title(nams{tsub});
                            axis xy square off
                        end
                        print(tempH,'-dpsc2',tempoutfile);
                        delete(tempH);
                    end

                    if strcmp(txt,sprintf('Superimposed on %s K-weighted contrast',Ktype));
                        fprintf('\nSaving K-weighted contrast')
                        V=spm_vol(fullfile(indir,effects{e},'spmF_0001.img'));
                        V.fname=fullfile(indir,effects{e},[Ktype 'K.img']);
                        spm_write_vol(V,Y);
                    end

                else
                    % factor just has two levels so beta file shows the contrast of
                    % the levels
                    if ~isempty(regexp(timewin,'saverage','ONCE'));continue; end
                    try
                        cfile=fullfile(indir,effects{e},'beta_0001.img');
                        Y=spm_read_vols(spm_vol(cfile));
                        txt='Superimposed on contrast';
                    catch error('Failed to find %s',cfile)
                    end
                end


                if ~isempty(regexp(timewin,'saverage','ONCE'));continue; end

                imagesc(Y'); hold on
                % add significance contours
                H=[];
                if exist(fullfile(indir,effects{e},'0p05_none.img'),'file')
                    Sunc=spm_read_vols(spm_vol(fullfile(indir,effects{e},'0p05_none.img')));
                    if any(Sunc(:))
                        [C H(1)]=contour(Sunc',[0 0]); set(H(1),'linecolor','k','linewidth',2,'linestyle',':');
                    else txt=['(NS uncorrected) ' txt];
                    end
                end
                if exist(fullfile(indir,effects{e},'0p05_FDR.img'),'file')
                    Sfdr=spm_read_vols(spm_vol(fullfile(indir,effects{e},'0p05_FDR.img')));
                    [C H(2)]=contour(Sfdr',[0 0]); set(H(2),'linecolor','k','linewidth',1);
                end
                if exist(fullfile(indir,effects{e},'0p05_FWE.img'),'file')
                    Sfdr=spm_read_vols(spm_vol(fullfile(indir,effects{e},'0p05_FWE.img')));
                    [C H(3)]=contour(Sfdr',[0 0]); set(H(3),'linecolor','k','linewidth',3);
                end
                % would like these in white, but sometimes they are turned black
                % when printing so then don't corespond to the legend which stays
                % white!

                %colormap jet
                %load /imaging/dm01/MEG/aaMEG/redblue2.mat; colormap(map)
                colormap(flipud(lbmap(128,'RedBlue')))
                % add colour scale,
                lim=max(abs(caxis));
                caxis([-lim lim]); %symmetrical about zero
                pos{ct}=get(sp{ct},'position');
                cb{ct}=colorbar;

                if ~tf
                    axis xy image
                    set(gca,'xTick',[],'yTick',[])
                    set(sp{ct},'position',pos{ct}); drawnow
                    cbpos=get(cb{ct},'position');
                    set(cb{ct},'position',[pos{ct}(1)+pos{ct}(3) cbpos(2) cbpos(3)/2 cbpos(4)])
                else
                    axis xy square
                    wt=regexp(indir,'ST?t\d_([^_]*)_[^_]*$','tokens');
                    ch=spm_select('FPList',aas_getsesspath(aap,1,1),[wt{1}{1} '.*mat']);
                    load(deblank(ch(1,:)));
                    xlims=xlim;
                    step=fix(range(xlims)/8);
                    ticks=xlims(1):step:xlims(2);
                    labels=((-D.events.start:step:D.events.stop)+1)/D.Radc*1000; % ms from trigger
                    labels=sort(unique([labels D.events.stop/D.Radc*1000]));
                    set(sp{ct},'xtick',ticks,'xticklabel',labels );
                    hold on
                    Z=plot(repmat(D.events.start,[1 2]),ylim,'k--');
                    if ct==1; ylabel('Hz'); end
                    cclim=caxis;
                    clim=[min(cclim(1),clim(1)) max(cclim(2),clim(2))];
                end

                load(fullfile(indir,effects{e},'SPM.mat'));
                title(sprintf('%s (N=%g)',chantype,SPM.nscan/SPM.factor.levels))

                if tf; continue; end

                % add sensors
                switch chantype
                    case {'eeg'}
                        load /imaging/local/spm/spm5/EEGtemplates/draft_eeg_montage_djm_mean.mat
                        topchan=size(Cpos,2);
                    case {'mags','grds','lats','longs'}
                        load /imaging/local/spm/spm5/EEGtemplates/FIF306_setup.mat
                        topchan=102;
                end
                trueCpos=Cpos*size(Y,1);
                sens=plot(trueCpos(1,1:topchan),trueCpos(2,1:topchan),'.');
                set(sens,'color',[.5 .5 .5])

                if ct==length(chantypes)
                    % add legend
                    leglabs={'p<0.05 uncorrected','p<0.05 FDR-corrected','p<0.05 FWE-corrected'};
                    leg=legend(H(H>0),leglabs(H>0),'Location','southOutside');
                    legend boxoff
                end

                switch chantype
                    case {'mags'}
                        set(get(cb{ct},'yLabel'),'string','fT')
                    case {'eeg'}
                        set(get(cb{ct},'yLabel'),'string','uV')
                    otherwise
                        set(get(cb{ct},'yLabel'),'string','a.u.')
                end

            end % next chantype

            if ~isempty(regexp(timewin,'saverage','ONCE'));continue; end

            if tf % equate colorbars?
                % can't get it to work
                for ct=1:length(chantypes)
                    % axis(sp{ct}); caxis(clim);
                    set(sp{ct},'position',pos{ct}); drawnow
                    cbpos=get(cb{ct},'position');
                    set(cb{ct},'position',[pos{ct}(1)+pos{ct}(3) cbpos(2) cbpos(3)/2 cbpos(4)])
                    f=str2num(aap.tasksettings.aamod_emeg_timefrequency.frequencies);
                    set(sp{ct},'ytickLabel',f(get(sp{ct},'yTick')));
                end
            end

            % add title
            annotation(F,'textbox',[0 0.88 1 .12],'String', ...
                {sprintf('Contours show significant F test for %s [%s]',effects{e},timewin),txt}, ...
                'edge','none');
            % print
            spm_print(summaryfile);
            print(F,'-dpng',outfile);
            clf(F);
        end % already done?
    end % next effect
end % next directory
fprintf('\ndone.\n')
