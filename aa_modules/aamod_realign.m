% AA module - realignment
% [aap,resp]=aamod_realign(aap,task,subj)
% Motion correction of EPI images (a.k.a. realignment) using SPM5
% Rhodri Cusack MRC CBU 2004-6 based on original by Matthew Brett
% Modified by Rik Henson 2006-8 to accept reslice "which" option
% 	(plus more defaults can be passed)
% Tibor Auer MRC CBU Cambridge 2012-2013

function [aap,resp]=aamod_realign(aap,task,subj)

resp='';

switch task
    case 'report' % reformatted [TA]
        mvmean=[];
        mvmax=[];
        mvstd=[];
        mvall=[];
        nsess=length(aap.acq_details.sessions);
        
        qq=[];
        
        aap = aas_report_add(aap,subj,'<table><tr>');
        for sess=1:nsess
            aap = aas_report_add(aap,subj,'<td>');
            aap = aas_report_add(aap,subj,['<h3>Session: ' aap.acq_details.sessions(sess).name '</h3>']);
            fn = fullfile(aas_getsubjpath(aap,subj),['diagnostic_aamod_realign_' aap.acq_details.sessions(sess).name '.jpg']);
            
            % Custom plotting [TA]
            mv = aas_plot_realignPars(aap,subj,sess,~exist(fn,'file'));
            
            aap.report.mvmax(subj,sess,:)=max(mv);
            % mvmean(sess,:)=mean(mv);
            % mvstd(sess,:)=std(mv);
            % mvall=[mvall;mv];
            aap=aas_report_addimage(aap,subj,fn);
            
            aap = aas_report_add(aap,subj,'<h4>Movement maximums</h4>');
            aap = aas_report_add(aap,subj,'<table cellspacing="10">');
            aap = aas_report_add(aap,subj,sprintf('<tr><td align="right">Sess</td><td align="right">x</td><td align="right">y</td><td align="right">z</td><td align="right">rotx</td><td align="right">roty</td><td align="right">rotz</td></tr>',sess));
            aap = aas_report_add(aap,subj,sprintf('<tr><td align="right">%d</td>',sess));
            aap = aas_report_add(aap,subj,sprintf('<td align="right">%8.3f</td>',aap.report.mvmax(subj,sess,:)));
            aap = aas_report_add(aap,subj,sprintf('</tr>',sess));
            aap = aas_report_add(aap,subj,'</table>');
            
            aap = aas_report_add(aap,subj,'</td>');
        end;
        aap = aas_report_add(aap,subj,'</tr></table>');
        
        varcomp=mean((std(mvall).^2)./(mean(mvstd.^2)));
        aap = aas_report_add(aap,subj,'<h3>All variance vs. within session variance</h3><table><tr>');
        aap = aas_report_add(aap,subj,sprintf('<td>%8.3f</td>',varcomp));
        aap = aas_report_add(aap,subj,'</tr></table>');
        
        aap=aas_report_addimage(aap,subj,fullfile(aas_getsubjpath(aap,subj),'diagnostic_aamod_realign.jpg'));

		% Summary in case of more subjects [TA]
        if (subj > 1) && (subj == numel(aap.acq_details.subjects)) % last subject
            
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
        for sess = aap.acq_details.selected_sessions %
            % get files from stream
            imgs(sess) = {aas_getimages_bystream(aap,subj,sess,'epi');};
        end
        
        % [AVG] This will ensure that any printing commands of SPM are done in the subject directory...
        cd(aas_getsubjpath(aap,subj))
        
        % Run the realignment
        spm_realign(imgs);
        
        if (~isdeployed)
            % Save graphical output
            try figure(spm_figure('FindWin', 'Graphics')); catch; figure(1); end
            print('-djpeg','-r150',fullfile(aas_getsubjpath(aap,subj),'diagnostic_aamod_realign'));
        end
        
        % Sessionwise custom plot [TA]
        for sess = aap.acq_details.selected_sessions
            aas_plot_realignPars(aap,subj,sess,true);
        end
        
        % Run the reslicing
        spm_reslice(imgs, resFlags);
        
        %% Describe outputs
        movPars = {};
        for sess = aap.acq_details.selected_sessions
            aas_log(aap,0,sprintf('Working with session %d', sess))
            
            rimgs=[];
            for k=1:size(imgs{sess},1);
                [pth nme ext]=fileparts(imgs{sess}(k,:));
                
                % If we don't reslice the images after realignment, don't
                % change the prefix of the images in our output stream
                %   - cwild 03/06/12
                if aap.tasklist.currenttask.settings.reslicewhich == 0
                    rimgs=strvcat(rimgs,fullfile(pth,[nme ext]));
                else
                    rimgs=strvcat(rimgs,fullfile(pth,['r' nme ext]));
                end
            end
            sessdir = aas_getsesspath(aap,subj,sess);
            aap = aas_desc_outputs(aap,subj,sess,'epi',rimgs);
            
            % Get the realignment parameters...
            fn=dir(fullfile(pth,'rp_*.txt'));
            outpars = fullfile(pth,fn(1).name);
            % Add it to the movement pars...
            movPars = [movPars outpars];
            
            aap = aas_desc_outputs(aap,subj,sess,'realignment_parameter', outpars);
            
            if find(sess==aap.acq_details.selected_sessions) == 1 % [AVG!]
                % mean only for first session
                fn=dir(fullfile(pth,'mean*.nii'));
                aap = aas_desc_outputs(aap,subj,'meanepi',fullfile(pth,fn(1).name));
            end
        end
        
        %% DIAGNOSTICS
        mriname = aas_prepare_diagnostic(aap,subj);
        
        aas_realign_graph(movPars)
        print('-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' mriname '_MP.jpeg']));
        
        cd(startingDir);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
end
