% AA module - tsdiffana - tool to assess time series variance
% Rhodri Cusack MRC CBU Cambridge Aug 2004
% subj=subject num
% sess=session num
% Original code written by Doris Eckstein and James Rowe
% Improved (and hopefully not broken) by Rhodri Cusack and Karolina
% Moutsopoulou Jun 2008

function [aap,resp]=aamod_listspikes(aap,task,subj,sess)

resp='';

switch task
    case 'report'
        spl = load(aas_getfiles_bystream(aap,subj,sess,'listspikes'));
        aap.report.html=strcat(aap.report.html,sprintf('Spikes %d   Moves %d<br>',size(spl.TSspikes,1),size(spl.Mspikes,1)));
    case 'doit'
        
        subjname = aas_prepare_diagnostic(aap,subj);
        
        try close(2); catch; end
        figure(2)
        
        % lists all scans with high deviation from mean, based on timediff.mat file
        % created in tsdiffana.m, and on rp*.txt file, created by
        % spm_realign.
        
        % Load the movement parameters
        Mfn = cellstr(aas_getfiles_bystream(aap,subj,sess,'realignment_parameter')); % aas_movPars(aap,subj, [1 0 0; 0 0 0]);
        rp = spm_load(Mfn{strcmp(spm_file(Mfn, 'ext'),'txt')});
        nsess = length(aap.acq_details.sessions);
        
        
        % Load up differnces through time as produced by tsdiffana
        tdfn = aas_getimages_bystream(aap,subj,sess,'tsdiffana');
        
        try
            qa = load (tdfn);
        catch
            aas_log(aap,1,sprintf('%s not found: Please run tsdiffana first',tdfn));
        end
        
        xyzlimit = aap.tasklist.currenttask.settings.xyzlimit;
        rotlimit_radians=aap.tasklist.currenttask.settings.rotlimit_degrees*pi/180;
        
        %% Now find big changes from one image to the next
        %  qa.qa.global.diff = mean (across voxels) of square difference between one volume and the next
        %  qa.qa.global.mean = mean global value across an image
        
        tm = qa.qa.global.diff/(mean(qa.qa.global.mean).^2); % RC/KM added .^2 16/6/2008
        
        switch aap.tasklist.currenttask.settings.tmbaseline
            case 'zero'
                Btm = 0;
            case 'mean'
                Btm = mean(tm);
            case 'median'
                Btm = median(tm);
            case 'smooth'
                Btm = smooth(1:length(tm),tm,0.1,'rloess');
        end
        
        % Residuals of the line
        Rtm = tm - Btm;
        
        switch aap.tasklist.currenttask.settings.tmmode
            case 'absolute'
                tmlimit = aap.tasklist.currenttask.settings.tmlimit;
            case 'std'
                tmlimit = aap.tasklist.currenttask.settings.tmlimit * std(Rtm);
            case 'rstd'
                % Robust std using only bottom 99% of distribution
                oRtm = sort(Rtm);
                oRtm = oRtm(1:round(0.99*length(oRtm)));
                tmlimit = aap.tasklist.currenttask.settings.tmlimit * std(oRtm);
        end
        
        badimages=[false; (Rtm > tmlimit)];
        
        TSspikes=[find(badimages),tm(badimages(2:end)),qa.qa.slice.diff(badimages(2:end))];
        
        %% Now find big movements
        % shift to sync with scan number
        rpdiff = [zeros(1,6); diff(rp)];
        absrpdiff = abs(rpdiff);
        
        badTspikes = any(absrpdiff(:,1:3) > xyzlimit,2);
        badRspikes = any(absrpdiff(:,4:6) > rotlimit_radians,2);
        badMspikes= badTspikes | badRspikes;
        
        Mspikes=[find(badMspikes), rpdiff(badMspikes,:)];
        
        %% DIAGNOSTIC
        
        subplot(nsess,4, (sess - 1) * 4 + 1)
        hold off
        plot(tm, 'b.')
        hold on
        plot(TSspikes(:,1), TSspikes(:,2), 'ko')
        title(sprintf('Sess %d \t Spikes: %d\n', sess, size(TSspikes,1)))
        
        subplot(nsess,4, (sess - 1) * 4 + 2)
        hist(Rtm,50)
        title(sprintf('Distribution of the tm data underlying spikes'))
        
        subplot(nsess,4, (sess - 1) * 4 + 3)
        hold off
        plot(rpdiff(:,1:3))
        hold on
        plot(find(badTspikes), rpdiff(badTspikes,1:3), 'ko')
        title(sprintf('Sess %d \t Translations: %d\n', sess, sum(badTspikes,1)))
        
        subplot(nsess,4, (sess - 1) * 4 + 4)
        hold off
        plot(rpdiff(:,4:6))
        hold on
        plot(find(badRspikes), rpdiff(badRspikes,4:6), 'ko')
        title(sprintf('Sess %d \t Rotations: %d\n', sess, sum(badRspikes,1)))
        
        
        %% Save things
        aas_log(aap,false,sprintf('Sess %d \t Spikes: %d; Moves: %d', sess, size(TSspikes,1), size(Mspikes,1)))
        
        SPfn = fullfile(aas_getsesspath(aap,subj,sess),sprintf('spikesandMspikes.mat'));
        save(SPfn, 'TSspikes', 'Mspikes');
        
        % Save the time differences
        aap = aas_desc_outputs(aap,subj,sess, 'listspikes', SPfn);
        
        
        %% Save graphical output to common diagnostics directory
        set(2,'Renderer','zbuffer');
        print(2,'-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' subjname '.jpeg']));
        try close(2); catch; end
    case 'checkrequirements'

    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
