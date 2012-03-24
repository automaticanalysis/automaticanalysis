% AA module - tsdiffana - tool to assess time series variance
% Rhodri Cusack MRC CBU Cambridge Aug 2004
% i=subject num
% j=session num
% Original code written by Doris Eckstein and James Rowe
% Improved (and hopefully not broken) by Rhodri Cusack and Karolina
% Moutsopoulou Jun 2008

function [aap,resp]=aamod_listspikes(aap,task,i,j)

resp='';

switch task
    case 'domain'
        resp='session';   % this module needs to be run once per session
    case 'description'
        resp='Run list spikes';
    case 'summary'
        resp='List spikes\n';
    case 'report'
        dirn = aas_getsesspath(aap,i,j);
        spfn = fullfile(dirn,'spikesandmoves.mat');
        load(spfn);
        aap.report.html=strcat(aap.report.html,sprintf('Spikes %d   Moves %d<br>',size(spikes,1),size(moves,1)));
    case 'doit'

        % lists all scans with high deviation from mean, based on timediff.mat file
        % created in tsdiffana.m, and on rp*.txt file, created by spm_realign.
        % i - subject num
        % j - session num
        
        % change upper limits here
        tmlimit = 0.05;            % threshold for mean square image difference, scaled by globals ^ 2
        xyzlimit = .3;          % translation limit from image to image in mm
        rotlimit_degrees = 2; % rotation limit in degrees

%         % These are the default values
%         tmlimit=aap.spmanalysis.spikesandmoves.spikethreshold; % difference between images (proportion of (sum squared diff)/(globals squared))
%         xyzlimit = aap.spmanalysis.spikesandmoves.translationthreshold; %mm
%         rotlimit_degrees = aap.spmanalysis.spikesandmoves.rotationthreshold; % degrees

        % Load up differnces through time as produced by tsdiffana
        dirn = aas_getsesspath(aap,i,j);
        tdfn = fullfile(dirn,'timediff.mat');
        try
            load (tdfn, 'td', 'globals', 'slicediff');
        catch
            aas_log(aap,1,sprintf('%s not found: Please run tsdiffana first',tdfn));
        end

        % Now find big changes from one image to the next
        %  td= mean (across voxels) of square difference between one volume and the next
        %  globals= mean global value across an image

        tm = td/(mean(globals).^2); % RC/KM added .^2 16/6/2008
        badimages=[false; (tm > tmlimit)] ;
        spikes=[find(badimages),tm(badimages(2:end)),slicediff(badimages(2:end))];


        % Now find big movements
        rotlimit_radians=rotlimit_degrees*pi/180;
        rpfn = spm_select('List',dirn,'^rp.*txt');
        moves = [];
        if (numel(rpfn) == 0),
            aas_log(aap,0,sprintf('^rp.*txt not found: No movement analysis done'));
        else
            rp=spm_load(fullfile(dirn,rpfn));
            % shift to sync with scan number
            rpdiff = [zeros(1,6); diff(rp)];
            absrpdiff = abs(rpdiff);

            badmoves=(any(absrpdiff(:,1:3) > xyzlimit,2)) | (any(absrpdiff(:,4:6) > rotlimit_radians,2));
            moves=[find(badmoves),rpdiff(badmoves,:)];
        end;



        if (isfield(aap.tasklist.currenttask.extraparameters,'filesuffix'))
            filesuffix=aap.tasklist.currenttask.extraparameters.filesuffix;
        else
            filesuffix=[];
        end;
        spfn = fullfile(dirn,sprintf('spikesandmoves%s.mat',filesuffix));
        save(spfn, 'spikes', 'moves');

    case 'checkrequirements'

    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
