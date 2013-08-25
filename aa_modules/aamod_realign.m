% AA module - realignment
% [aap,resp]=aamod_realign(aap,task,i)
% Motion correction of EPI images (a.k.a. realignment) using SPM5
% Rhodri Cusack MRC CBU 2004-6 based on original by Matthew Brett
% Modified by Rik Henson 2006-8 to accept reslice "which" option
% 	(plus more defaults can be passed)
% Tibor Auer MRC CBU Cambridge 2012-2013

function [aap,resp]=aamod_realign(aap,task,i)

resp='';

switch task
    case 'report' % reformatted [TA]
        mvmean=[];
        mvmax=[];
        mvstd=[];
        mvall=[];
        nsess=length(aap.acq_details.sessions);
        
        qq=[];
		
		aap = aas_report_add(aap,i,'<table><tr>');
        for sess=1:nsess
            aap = aas_report_add(aap,i,'<td>');
            aap = aas_report_add(aap,i,['<h3>Session: ' aap.acq_details.sessions(sess).name '</h3>']);
            fn = fullfile(aas_getsubjpath(aap,i),['diagnostic_aamod_realign_' aap.acq_details.sessions(sess).name '.jpg']);
           
			% Custom plotting [TA]
            mv = plot_parameters(aap,i,sess,~exist(fn,'file'));

            aap.report.mvmax(i,sess,:)=max(mv);           
            % mvmean(sess,:)=mean(mv);
            % mvstd(sess,:)=std(mv);
            % mvall=[mvall;mv];            
            aap=aas_report_addimage(aap,i,fn);
            
            aap = aas_report_add(aap,i,'<h4>Movement maximums</h4>');
            aap = aas_report_add(aap,i,'<table cellspacing="10">');
            aap = aas_report_add(aap,i,sprintf('<tr><td align="right">Sess</td><td align="right">x</td><td align="right">y</td><td align="right">z</td><td align="right">rotx</td><td align="right">roty</td><td align="right">rotz</td></tr>',sess));
            aap = aas_report_add(aap,i,sprintf('<tr><td align="right">%d</td>',sess));
            aap = aas_report_add(aap,i,sprintf('<td align="right">%8.3f</td>',aap.report.mvmax(i,sess,:)));
            aap = aas_report_add(aap,i,sprintf('</tr>',sess));
            aap = aas_report_add(aap,i,'</table>');
            
            aap = aas_report_add(aap,i,'</td>');
        end;
        aap = aas_report_add(aap,i,'</tr></table>');
        
        varcomp=mean((std(mvall).^2)./(mean(mvstd.^2)));
        aap = aas_report_add(aap,i,'<h3>All variance vs. within session variance</h3><table><tr>');
        aap = aas_report_add(aap,i,sprintf('<td>%8.3f</td>',varcomp));
        aap = aas_report_add(aap,i,'</tr></table>');
        
        aap=aas_report_addimage(aap,i,fullfile(aas_getsubjpath(aap,i),'diagnostic_aamod_realign.jpg'));

		% Summary in case of more subjects [TA]
        if (i > 1) && (i == numel(aap.acq_details.subjects)) % last subject
            meas = {'Trans - x','Trans - y','Trans - z','Pitch','Roll','Yaw'};
            for sess=1:nsess
                mvmax = squeeze(aap.report.mvmax(:,sess,:));
                boxplot(mvmax,'label',meas);
                boxValPlot = getappdata(getappdata(gca,'boxplothandle'),'boxvalplot');
                fn = fullfile(aas_getstudypath(aap),['diagnostic_aamod_realignunwarp_' aap.acq_details.sessions(sess).name '.jpg']);
                print('-djpeg','-r75',fn);
                close(gcf);
                
                aap = aas_report_add(aap,'moco','<td>');
                aap = aas_report_add(aap,'moco',['<h3>Session: ' aap.acq_details.sessions(sess).name '</h3>']);
                aap=aas_report_addimage(aap,'moco',fn);
                
                for ibp = 1:numel(meas)
                    bp = boxValPlot(ibp,:);
                    subjs = ' None';
                    if bp.numFiniteHiOutliers                        
                        subjs = [' ' num2str(sort(cell2mat(bp.outlierrows)'))];
                    end
                    aap = aas_report_add(aap,'moco',sprintf('<h4>Outlier(s) in %s:%s</h4>',meas{ibp},subjs));
                end
                
                aap = aas_report_add(aap,'moco','</td>');
            end
        end        
    case 'doit'
        % Get realignment defaults
        defs = aap.spm.defaults.realign;
       
	   
        % Note starting directory so we can get back here in the end
        startingDir = pwd;
        
		% Flags to pass to routine to calculate realignment parameters
        % (spm_realign)
        reaFlags = struct(...
            'quality', defs.estimate.quality,...  % estimation quality
            'fwhm', defs.estimate.fwhm,...        % smooth before calculation
            'rtm', defs.estimate.rtm,...          % whether to realign to mean
            'interp', defs.estimate.interp,...    % interpolation type
            'wrap', defs.estimate.wrap,...        % wrap in (x) y (z)
            'sep', defs.estimate.sep...          % interpolation size (separation)
            );
        
        % Flags to pass to routine to create resliced images
        % (spm_reslice)
        resFlags = struct(...
            'interp', defs.write.interp,...       % interpolation type
            'wrap', defs.write.wrap,...           % wrapping info (ignore...)
            'mask', defs.write.mask,...           % masking (see spm_reslice)
            'which', aap.tasklist.currenttask.settings.reslicewhich,...     % what images to reslice
            'mean', aap.tasklist.currenttask.settings.writemean);           % write mean image
        
        clear imgs;
        for j = aap.acq_details.selected_sessions %
            % get files from stream
            imgs(j) = {aas_getimages_bystream(aap,i,j,'epi');};
        end
        
        % [AVG] This will ensure that any printing commands of SPM are done in the subject directory...
        cd(aas_getsubjpath(aap,i))
        
        % Run the realignment
        spm_realign(imgs);
        
        if (~isdeployed)
            % Save graphical output
            try figure(spm_figure('FindWin', 'Graphics')); catch; figure(1); end;
            print('-djpeg','-r75',fullfile(aas_getsubjpath(aap,i),'diagnostic_aamod_realign'));
        end;
		
		% Sessionwise custom plot [TA]
        for sess = aap.acq_details.selected_sessions
            plot_parameters(aap,i,sess,true);
        end

		
        % Run the reslicing
        spm_reslice(imgs, resFlags);
        
        % Describe outputs
        for j = aap.acq_details.selected_sessions
            rimgs=[];
            for k=1:size(imgs{j},1);
                [pth nme ext]=fileparts(imgs{j}(k,:));
                
                % If we don't reslice the images after realignment, don't
                % change the prefix of the images in our output stream
                %   - cwild 03/06/12
                if aap.tasklist.currenttask.settings.reslicewhich == 0            
                    rimgs=strvcat(rimgs,fullfile(pth,[nme ext]));
                else
                    rimgs=strvcat(rimgs,fullfile(pth,['r' nme ext]))
                end
            end;
            sessdir=aas_getsesspath(aap,i,j);
            aap = aas_desc_outputs(aap,i,j,'epi',rimgs);
            
            fn=dir(fullfile(pth,'rp_*.txt'));
            aap = aas_desc_outputs(aap,i,j,'realignment_parameter',fullfile(pth,fn(1).name));
            
            if (j==1)
                % mean only for first session
                fn=dir(fullfile(pth,'mean*.nii'));
                aap = aas_desc_outputs(aap,i,1,'meanepi',fullfile(pth,fn(1).name));
            end;
            
        end;
        
        cd(startingDir);
		
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
end

function RealPar = plot_parameters(aap,subj,sess,toDisp)
% Based on SPM Realign

% Threshold for excessive movement
QA_TRANSL = 2;
QA_ROT = 8;

P = spm_vol(aas_getimages_bystream(aap,subj,sess,'epi'));
if length(P)<2, return; end;
Params = zeros(numel(P),12);
for i=1:numel(P),
    Params(i,:) = spm_imatrix(P(i).mat/P(1).mat);
end
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
end