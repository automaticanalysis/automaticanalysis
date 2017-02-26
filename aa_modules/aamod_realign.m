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
        if subj == 1 % init summary
            aap.report.(mfilename).selected_sessions = zeros(1,0);
            aap.report.(mfilename).mvmax = nan(aas_getN_bydomain(aap,'subject'),aas_getN_bydomain(aap,'session',1),6);
        end
        aap.report.(mfilename).selected_sessions = union(aap.report.(mfilename).selected_sessions,aap.acq_details.selected_sessions);
        
        mvmean=[];
        mvmax=[];
        mvstd=[];
        mvall=[];
        nsess=length(aap.acq_details.selected_sessions);
        
        qq=[];
        
        aap = aas_report_add(aap,subj,'<table><tr>');
        for sess=aap.acq_details.selected_sessions
%             if sess > aas_getN_bydomain(aap,'session',subj), break; end
            aap = aas_report_add(aap,subj,'<td>');
            aap = aas_report_add(aap,subj,['<h3>Session: ' aap.acq_details.sessions(sess).name '</h3>']);
            fn = fullfile(aas_getsubjpath(aap,subj),['diagnostic_aamod_realign_' aap.acq_details.sessions(sess).name '.jpg']);
            
            par = cellstr(aas_getfiles_bystream(aap,subj,sess,'realignment_parameter'));
            parind = cell_index(par,'.txt');
            mv = load(par{parind});

            if ~exist(fn,'file')
                if isfield(aap.tasklist.currenttask.settings,'mfp') && aap.tasklist.currenttask.settings.mfp.run
                    mw_mfp_show(aas_getsesspath(aap,subj,sess));
                    movefile(...
                        fullfile(aas_getsesspath(aap,subj,sess),'mw_motion.jpg'),fn);
                else
                    f = aas_realign_graph(par{parind});
                    print('-djpeg','-r150','-noui',...
                        fullfile(aas_getsubjpath(aap,subj),...
                        ['diagnostic_aamod_realignunwarp_' aap.acq_details.sessions(sess).name '.jpg'])...
                        );
                    close(f);
                end
            end

            aap.report.(mfilename).mvmax(subj,sess,:)=max(mv);
            % mvmean(sess,:)=mean(mv);
            % mvstd(sess,:)=std(mv);
            % mvall=[mvall;mv];
            aap=aas_report_addimage(aap,subj,fn);
            
            aap = aas_report_add(aap,subj,'<h4>Movement maximums</h4>');
            aap = aas_report_add(aap,subj,'<table cellspacing="10">');
            aap = aas_report_add(aap,subj,sprintf('<tr><td align="right">Sess</td><td align="right">x</td><td align="right">y</td><td align="right">z</td><td align="right">rotx</td><td align="right">roty</td><td align="right">rotz</td></tr>',sess));
            aap = aas_report_add(aap,subj,sprintf('<tr><td align="right">%d</td>',sess));
            aap = aas_report_add(aap,subj,sprintf('<td align="right">%8.3f</td>',aap.report.(mfilename).mvmax(subj,sess,:)));
            aap = aas_report_add(aap,subj,sprintf('</tr>',sess));
            aap = aas_report_add(aap,subj,'</table>');
            
            aap = aas_report_add(aap,subj,'</td>');
        end;
        aap = aas_report_add(aap,subj,'</tr></table>');
        
        varcomp=mean((std(mvall).^2)./(mean(mvstd.^2)));
        aap = aas_report_add(aap,subj,'<h3>All variance vs. within session variance</h3><table><tr>');
        aap = aas_report_add(aap,subj,sprintf('<td>%8.3f</td>',varcomp));
        aap = aas_report_add(aap,subj,'</tr></table>');
        
		% Summary in case of more subjects [TA]
        if (subj > 1) && (subj == numel(aap.acq_details.subjects)) % last subject            
            meas = {'Trans - x','Trans - y','Trans - z','Pitch','Roll','Yaw'};
            for sess=aap.report.(mfilename).selected_sessions
				fn = fullfile(aas_getstudypath(aap),['diagnostic_aamod_realign_' aap.acq_details.sessions(sess).name '.jpg']);
                
                mvmax = squeeze(aap.report.(mfilename).mvmax(:,sess,:));
                f = figure; boxplot(mvmax,'label',meas);
                boxValPlot = getappdata(getappdata(gca,'boxplothandle'),'boxvalplot');
                set(f,'Renderer','zbuffer');
                if ~exist(fn,'file'), print(f,'-djpeg','-r150',fn); end
                close(f);
                
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
        elseif numel(aap.acq_details.subjects) == 1
            aap = aas_report_add(aap,'moco','<h4>No summary is generated: there is only one subject in the pipeline</h4>');
        end
    case 'doit'
        % Get realignment defaults from the XML!
        reaFlags = aap.tasklist.currenttask.settings.eoptions;
        if strcmp(reaFlags.weight,'''''') % empty
            reaFlags.weight = {};
        else
            reaFlags.weight = {reaFlags.weight};
        end        
        
        resFlags = aap.tasklist.currenttask.settings.roptions;
        resFlags.which = [aap.tasklist.currenttask.settings.reslicewhich ...
            aap.tasklist.currenttask.settings.writemean];     
        
        % Note starting directory so we can get back here in the end
        startingDir = pwd;
        
        imgs = {};
        for sess = aap.acq_details.selected_sessions
            % get files from stream
            imgs{end+1} = aas_getfiles_bystream(aap, subj, sess, 'epi');
        end
        
        % Check if we are using a weighting image
        if any(strcmp(aap.tasklist.currenttask.inputstreams.stream, 'weightingImage'))
            
            wImgFile = aas_getfiles_bystream(aap,subj,'weightingImage');
            wVol = spm_vol(wImgFile);

            aas_log(aap,false,sprintf('Realignment is going to be weighted with: %s', wVol.fname));
            
            % Use the first EPI as a space reference 
            rVol = spm_vol(imgs{1}(1,:));
            
            % Check if the dimensions and the orientation of the weighting
            % image match that of the first EPI in the data set.  If not,
            % we reslice the weighting image.
            if ( any(any(~(wVol.mat == rVol.mat))) ||  any(~(wVol.dim == rVol.dim)) )
                spm_reslice(strvcat(rVol.fname, wVol.fname), struct('which', 1, 'mean', 0, 'interp', 0, 'prefix', 'r'));
                [rPath rFile rExt] = fileparts(rVol.fname);
                [wPath wFile wExt] = fileparts(wImgFile);
                wImgFile = fullfile(wPath,['r' wFile wExt]);
                wVol = spm_vol(wImgFile);
            end
            
            if aap.tasklist.currenttask.settings.invertWeightingImage
                wY = spm_read_vols(wVol);
                wY = ~wY; % assuming binary weighting for now.
                spm_write_vol(wVol, wY);   
            end

            reaFlags.PW = wImgFile;     
        end       
        
        % [AVG] This will ensure that any printing commands of SPM are done in the subject directory...
        cd(aas_getsubjpath(aap,subj))
        
        % Run the realignment
        spm_realign(imgs, reaFlags);
        
        % Run the reslicing
        spm_reslice(imgs, resFlags);
        
        %% Describe outputs
        movPars = {};
        for sess = aap.acq_details.selected_sessions
            aas_log(aap,0,sprintf('Working with session %d: %s', sess, aap.acq_details.sessions(sess).name))
            
            rimgs=[];
            for k=1:size(imgs{aap.acq_details.selected_sessions==sess},1);
                [pth nme ext]=fileparts(imgs{aap.acq_details.selected_sessions==sess}(k,:));
                
                % If we don't reslice the images after realignment, don't
                % change the prefix of the images in our output stream
                %   - cwild 03/06/12
                if aap.tasklist.currenttask.settings.reslicewhich == 0
                    rimgs=strvcat(rimgs,fullfile(pth,[nme ext]));
                else
                    rimgs=strvcat(rimgs,fullfile(pth,['r' nme ext]));
                end
            end
            aap = aas_desc_outputs(aap,subj,sess,'epi',rimgs);
            
            % Get the realignment parameters...
            fn=dir(fullfile(pth,'rp_*.txt'));
            outpars = fullfile(pth,fn(1).name);
            
            % Sessionwise custom plot or MFP
            if isfield(aap.tasklist.currenttask.settings,'mfp') && aap.tasklist.currenttask.settings.mfp.run
                mw_mfp(outpars);
                fn=dir(fullfile(pth,'mw_mfp_*.txt'));
                outpars = fullfile(pth,fn(1).name);
                if strcmp(aap.options.wheretoprocess,'localsingle')
                    movefile(...
                        fullfile(aas_getsesspath(aap,subj,sess),'mw_motion.jpg'),...
                        fullfile(aas_getsubjpath(aap,subj),...
                        ['diagnostic_aamod_realign_' aap.acq_details.sessions(sess).name '.jpg'])...
                        );
                end
            else
                f = aas_realign_graph(outpars);
                print('-djpeg','-r150','-noui',...
                    fullfile(aas_getsubjpath(aap,subj),...
                    ['diagnostic_aamod_realign_' aap.acq_details.sessions(sess).name '.jpg'])...
                    );
                close(f);
            end
            
            aap = aas_desc_outputs(aap,subj,sess,'realignment_parameter', outpars);
            
            if sess==min(aap.acq_details.selected_sessions)
                % mean only for first session
                fn=dir(fullfile(pth,'mean*.nii'));
                aap = aas_desc_outputs(aap,subj,'meanepi',fullfile(pth,fn(1).name));
            end
            
            % Add it to the movement pars...
            movPars = [movPars outpars];
        end
        
        %% DIAGNOSTICS
        subjname = aas_prepare_diagnostic(aap,subj);
        
        f = aas_realign_graph(movPars);
        print('-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' subjname '_MP.jpeg']));
        close(f);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
end
