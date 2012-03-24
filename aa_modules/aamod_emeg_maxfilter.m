function [aap resp]=aamod_emeg_maxfilter(varargin)
% Runs MaxFilter from aap.
%
% When run from aa, format will be: aamod_maxfilter(aap,job,subject,block)
%
% Includes most options but not yet -headpos or advanced movecomp options.
% Saves data as 'short' 16 bit format suitable for SPM and BESA
% If MaxMove is used to align data to a reference fif, the HPI data
% is sometimes removed, so here it's also saved in a seperate file.
% If spm_eeg_rdata_FIF or spm_eeg_rdata_FIF_dm are updated, this script might
% need to be changed.
% Don't use -lpfilt as it introduces a spurious time shift.
% MaxFilter will crash if downsampling by a factor other than (1,) 2, 3 or 4
% Doesn't force
% The correct crosstalk and calibration files are important for
% maxfilter...is it possible to confirm these from the fif file?
%
% Danny Mitchell 03/03/08


%% check task settings, subject, block etc
[aap subblock doit resp settings]=aa_emeg_checktasksettings(mfilename('fullpath'),varargin);
if ~doit; return; end

%% look for raw data files from any date
clear files
subdir=fullfile(aap.directory_conventions.rawmegdatadir,aap.acq_details.subjects(subblock{1}).megname);
[junk dates]=spm_select('List', subdir, '');
if length(subblock)==2
    filt=regexprep(settings.InputFilter,'BLOCK',aap.acq_details.sessions(subblock{2}).name);
else filt=settings.InputFilter;
end
for d=3:size(dates,1)
    temp=spm_select('FPList',fullfile(subdir,deblank(dates(d,:))),filt);
    if ~isempty(temp);
        if exist('files','var'); files=cellstr(files,temp);
        else files=cellstr(temp);
        end
    end
end
try if exist(files{1},'file')~=2; files(1)=[]; end % could be'/'
catch files=[];
end
if isempty(files); aas_log(aap,1,sprintf('\nFound no data!\n Input filter is %s\n Looking in %s/* \n',settings.InputFilter,subdir)); end

%% move to output directory
blockdir=fullfile(aap.acq_details.root,aap.acq_details.subjects(subblock{1}).megname,aap.acq_details.sessions(subblock{2}).name);
cd(blockdir)

%% run task for each file
for f=1:length(files);
    [pth fn ext]=fileparts(files{f});

%% prepare output filenames, including prefix for SSS and MaxMove options
    if ~settings.UseSSS; prefix='NS_';
    else
        if settings.UseMaxST; prefix='ST_'; else prefix='S_'; end
        if settings.MaxMoveMethod<6;
            prefix=strrep(prefix,'_',sprintf('t%g_',settings.MaxMoveMethod));
        end
    end

    outfile=fullfile(blockdir,[prefix fn ext]);
    logfile=regexprep(outfile,'\.fif$','_MaxFilterLog.txt');
    logfile1=strrep(logfile,'.txt','1.txt');
   
    if exist(outfile,'file')==0 || settings.Overwrite==1;

        warning off all; try delete(outfile); catch end; warning on all % needed to overwrite if not forcing maxfilter

%% generate core maxfilter command
        if ~settings.UseSSS
            % just convert to short and downsample if requested
            cmd=sprintf('/neuro/bin/util/maxfilter -nosss -format short -v -f %s -o %s -ds %g', ...
                files{f},outfile, settings.DownSampleFactor);
            callmaxfilter(cmd,logfile1,aap)
        else
            % generate core MaxFilter execution command:
            cmd=sprintf('/neuro/bin/util/maxfilter -v -in %g -out %g -ctc %s -cal %s', ...
                settings.L_in, settings.L_out, ...
                settings.CrosstalkFile, settings.CalibrationFile);

            
%% backup digitised points in head coordinates (they're deleted 
% by -trans) and get digitised scalp points (not fids or EEG) in 
% head coordinates, removing outliers (>2sd from mean)
            HPI=strrep(outfile,'.fif','_HPI.mat');
            [co,ki, nu] = hpipoints(files{f});
            if exist(HPI,'file')==0 || settings.Overwrite==1;
                save(HPI, 'co', 'ki', 'nu');
            end
            headpoints=co(:,~mod(ki,2))';
            ulim=mean(headpoints)+2*std(headpoints);
            llim=mean(headpoints)-2*std(headpoints);
            xol=(headpoints(:,1)>ulim(1) | headpoints(:,1)<llim(1));
            yol=(headpoints(:,2)>ulim(2) | headpoints(:,2)<llim(2));
            zol=(headpoints(:,3)>ulim(3) | headpoints(:,3)<llim(3));
            headpoints(any([xol yol zol],2),:)=[];
            save(regexprep(outfile,'\.fif$','_headpoints.txt'),'-ASCII','headpoints')

%% get centre of best-fitting sphere to these points
            cmd_fit=['/neuro/bin/util/fit_sphere_to_points ' regexprep(outfile,'\.fif$','_headpoints.txt')];
            [status spherefit]=unix(cmd_fit);
            spherefit=str2num(spherefit)*1000; % vector to mm;
                    
%% generate MaxMove transformation command if requested
            % 1) None
            % 2) To first block of each block type per subject
            % 3) Neuromag default 1 step (sphere fit in head coordinates to device origin of 0,0,0)
            % 4) Neuromag default optimized (sphere fit in head coordinates to optimal device origin of 0,13,-6)
            % 5) Specify target fiff file (with expansion sphere at spherefit)
            clear co* ki nu
            if settings.MaxMoveMethod>1;
                if settings.MaxMoveMethod==2
                    MaxMove=sprintf('-trans %s -frame head -origin %g %g %g', files{1},spherefit(1),spherefit(2),spherefit(3));
                elseif settings.MaxMoveMethod==3
                    MaxMove=sprintf('-trans default -frame head -origin %g %g %g',spherefit(1),spherefit(2),spherefit(3)); % specify origin to avoid MaxFilter bug
                elseif settings.MaxMoveMethod==4
                    % Add transformation to maxfilter command, to align
                    % centre of head to optimal SSS origin for device
                    % (0,13, -6 in device coordinates). However, expansion 
                    % origin seems to be fixed to the device origin...
                    % If want to change this (e.g. keep it at head sphere
                    % origin), can do this in second step.
                    MaxMove=sprintf('-trans default -frame head -origin %g %g %g',spherefit(1),spherefit(2)-13,spherefit(3)+6);
                elseif settings.MaxMoveMethod==5 % to another fif file
                    MaxMove=sprintf('-trans %s -frame head -origin %g %g %g', settings.MaxMoveTarget,spherefit(1),spherefit(2),spherefit(3));
                end
            else MaxMove='';
            end
            
%% first pass: 
            % Use autobad to find static bad channels over 1st 20s (JT's
            % suggestion) - as can't use autobad with movecomp
            % Also create potential (temporary) target fif for MaxMove 
            if settings.AutoBad>0 || settings.MaxMoveMethod==4
                temptarget=regexprep(outfile,'.fif','_temptarget.fif');
                if exist(temptarget,'file') && settings.Overwrite; delete(temptarget); end
                if ~exist(temptarget,'file')
                    %cmd_1stpass=sprintf('%s -autobad on -skip 20 99999',cmd);
                    % All output local to each 1s buffer; occasional strange massive channel numbers
                    % Jason's text file collects each buffer in seperate rows

                    cmd_1stpass=sprintf('%s -f %s -o %s %s -autobad 20 -badlimit 7 -skip 20 99999',cmd,files{f},temptarget,MaxMove);
                    % Bad channels determined for initial chunk, then static bad
                    % channels defined and propogated to individual buffers;
                    % odd channel numbers don't appear
                    % Jason's text file collects the initially determined channels

                    fprintf('\nScanning first 20s for static bad channels, & creating -trans target. Please wait...\n')
                    callmaxfilter(cmd_1stpass,logfile1,aap)

%% send channel numbers from each buffer (per row) to a text file
                    badchanfile=strrep(logfile,'.txt','_badchans.txt');
                    if exist(badchanfile,'file'); delete(badchanfile); end
                    %getbad=sprintf('cat %s | sed -n ''/Detected/p'' | cut -f 5- -d '' '' > %s',logfile1,badchanfile);
                    getbad=sprintf('cat %s | sed -n ''/Static bad channels/p'' | cut -f 5- -d '' '' |uniq > %s',logfile1,badchanfile);
                    unix(getbad);
                    
                else fprintf('\n1st pass already completed.\n') 
                end
                
%% read consistently bad channels found by autobad
                fid = fopen(strrep(logfile,'.txt','_badchans.txt'),'rt');
                autobadchans = fgetl(fid);
                fclose(fid);
                if ischar(autobadchans); autobadchans=str2num(autobadchans);
                else autobadchans=[];
                end
                fprintf(' Found %g autobad channels (%s)',length(autobadchans),num2str(autobadchans))

%% add user specified bad chans from file if found
                acqfile=fullfile(aap.acq_details.root,'bad_chans.xls');
                if exist(acqfile,'file')
                    warning off all; [Numeric,Txt,Both]=xlsread(acqfile); warning on all
                    r=strmatch(aap.acq_details.subjects(subblock{1}).megname,Txt(:,1));
                    c=strmatch(regexprep(fn,'_raw',''),Txt(1,:));
                    val=Both{r,c};
                    if isnan(val); userbadchans=[];
                    elseif ischar(val); userbadchans=str2num(val);
                    else userbadchans=val;
                    end
                    fprintf('\n %g channels marked bad by user (%s)',length(userbadchans),num2str(userbadchans))
                else
                    fprintf('\nNo user-specified bad channel file found');
                    userbadchans=[];
                end
        
                % add to maxfilter command
                cmd=sprintf('%s -bad %s ',cmd,num2str(union(autobadchans,userbadchans)));
            end

%% add MaxST/autobad as requested, and MaxMove XOR movecomp
            if ~ischar(settings.MaxST_window) 
                % autobad not need ed with ST
                cmd_2ndpass=sprintf('%s -f %s -st %g -corr %g -autobad off', ...
                    cmd, files{f}, settings.MaxST_window, settings.MaxST_corr);
            else
                % just add autobad
                cmd_2ndpass=sprintf('%s -f %s -autobad %g',cmd,files{f},settings.AutoBad);
            end
            
            if settings.MaxMoveMethod==4
               % having created temptarget in the 'default' position but 
               % translated as desired, should now be able to use this as 
               % the target with the SSS origin at the origin of the head
               % sphere.
               MaxMove=sprintf('-trans %s -frame head -origin %g %g %g',temptarget,spherefit(1),spherefit(2),spherefit(3)); % spherefit(2)-15
            end
           
            if settings.UseMoveComp && settings.MaxMoveMethod
               % will do downsampling and MaxMove in a subsequent step... 
               % (specify SSS origin as spherefit) 
               movecomped=regexprep(outfile,'.fif','_movecomped.fif');
               if exist(movecomped,'file') && settings.Overwrite; delete(movecomped); end
               cmd_2ndpass=sprintf('%s -o %s -movecomp -autobad off -frame head -origin %g %g %g', ...
                   cmd_2ndpass,movecomped,spherefit(1),spherefit(2),spherefit(3));
            elseif settings.UseMoveComp
               % write to final destination file with movecomp
               % (specify SSS origin as spherefit) 
               cmd_2ndpass=sprintf('%s -o %s -ds %g -movecomp -autobad off -frame head -origin %g %g %g -format short', ...
                   cmd_2ndpass, outfile, settings.DownSampleFactor,spherefit(1),spherefit(2),spherefit(3));
            else
               % write to final destination file with trans 
               cmd_2ndpass=sprintf('%s -o %s -ds %g %s -format short',cmd_2ndpass, outfile, settings.DownSampleFactor,MaxMove);
            end
        
%% now call maxfilter again
        logfile2=strrep(logfile,'.txt','2.txt');
        if exist(outfile,'file'); delete(outfile); end 
        % ST is very slow: about 2 x acquisition time. (Although ST window length doesn't make huge difference.)
        % -autobad 1 is also fairly slow:
        % with neither ST nor autobad: about 0.25 x acquisition time 
        extratext=sprintf('%g user bad channels (%s) \n %g auto bad channels (%s) \n', ...
            length(userbadchans),num2str(userbadchans),length(autobadchans),num2str(autobadchans));
        if settings.UseMoveComp && exist(movecomped,'file')
            fprintf('\n2nd pass already completed.\n')
        else callmaxfilter(cmd_2ndpass,logfile2,aap,extratext);
        end
        
%% now call maxfilter for final time if wanted to trans after movecomp
        if settings.UseSSS && settings.UseMoveComp && settings.MaxMoveMethod
            cmd_3rdpass=sprintf('%s -f %s -o %s -ds %g %s -force -format short -autobad off ',cmd,movecomped,outfile,settings.DownSampleFactor,MaxMove);       
            logfile3=strrep(logfile,'.txt','3.txt');
            extratext=sprintf('%g user bad channels (%s) \n%g auto bad channels (%s) \n', ...
            length(userbadchans),num2str(userbadchans),length(autobadchans),num2str(autobadchans));
            callmaxfilter(cmd_3rdpass,logfile3,aap,extratext);
            delete(movecomped);
        end
        end % use SSS?
    end % already done?
    
    % work out size of rotation and translation applied
    if isnumeric(settings.MaxMoveMethod) && settings.MaxMoveMethod>1;
        T1=loadtrans(files{f}); % 'device to head1'
        T2=loadtrans(outfile); % 'device to head2'
        T12=T1\T2; % = inv(T1)*T2 = 'head1 to head2' ? % not sure why inv(T2)*T1 doesn't give opposite translations; could it be rounding errors??
        P=spm_imatrix(T12);
        transl=sqrt(sum(P(1:3).^2))*1000; % m to mm
        rot=acos((trace(T12)-2)/2) * 180/pi; % rad to deg; 
        % formula from Wikipedia:'rotation matrix', remembering that 
        % the 4x4 transformation matrix has an extra 1 on the diagonal 
        tits={'Transl','Rot'};
        vals={transl,rot};
        defs={'Total translation in mm','Total rotation in degrees'};
        aas_emeg_savestats(outfile,tits,vals,defs,settings.Overwrite);
    end

    fprintf('.');
end % next file

return

%%
function callmaxfilter(command,logfile,aap,extratext,showoutput)
if nargin<5; showoutput=true; end
if exist(logfile,'file'); delete(logfile); end;

diary(logfile);
if exist('extratext','var') && ~isempty(extratext); disp(extratext); end
disp(command);
tic;

if showoutput; status=unix(command);
else  [status output]=unix(command);
end

if status==6;
    % for some reason 1st tag may not be at time 0; if 1st
    % tag is at t>20s then all tags will be skipped and no
    % file will be written, so try without skipping.
    fprintf('\nTime of 1st tag may be > 20s; retrying without skipping any tags...');
    command=strrep(command,'-skip 20 99999','');
    delete(logfile);
    try
        status=unix(command);
    catch aas_log(aap,true,sprintf('\nFailed to run MaxFilter.\n'))
    end
elseif status>0

    diary off; 
    fid = fopen(logfile,'rt');
    tline=''; warnings='';
    while 1
        tline = fgetl(fid);
        if ~ischar(tline); break;
        elseif ~isempty(regexp(tline,'AUTOMATIC ANALYSIS FAILED','once')); break
        elseif ~isempty(regexp(tline,'Warning:','once'))
            warnings=tline;
        elseif ~isempty(warnings)
            warnings=char(warnings,tline);
        end
    end
    fclose(fid);

    if ~isempty(warnings)
        delete(logfile); diary(logfile);
        disp(warnings);
        fprintf('\nRetrying with force...\n');
        command=[command ' -force'];
        disp(command);
        status=unix(command);
    end
end

if status>0 % possibly no license for maxfilter??
    aas_log(aap,1,sprintf('\nFailed to run MaxFilter.\n'));
end

fprintf('\nProcessing took: %g minutes.\n',round(toc/60))
diary off

return