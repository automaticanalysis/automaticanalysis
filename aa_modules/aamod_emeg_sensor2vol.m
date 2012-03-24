function [aap resp]=aamod_emeg_sensor2vol(varargin)
% Danny Mitchell 02/04/08

%% check task settings, subject, block etc
[aap subblock doit resp settings]=aa_emeg_checktasksettings(mfilename('fullpath'),varargin);

%% CHECK REQUIREMENTS: load spreadsheet containing epoching instructions
epochfile=spm_select('FPList',fullfile(aas_getsesspath(aap,subblock{1},subblock{2}),'events'),'(E|e)pochs.*.xls');
if exist(epochfile,'file')~=2; epochfile=spm_select('FPList',fullfile(aas_getsubjpath(aap,subblock{1}),'events'),'(E|e)pochs.*\.xls'); end
if exist(epochfile,'file')~=2; epochfile=spm_select('FPList',fullfile(aap.acq_details.root,'events'),'(E|e)pochs.*.xls'); end
if exist(epochfile,'file')~=2; error('aa:EpochFileError', '\n Failed to find Epochs.*.xls, required for epoching continuous data. \n This should be placed in "events" folder at either study, subject, or session levels, with lower levels overriding higher ones. \n See examples/Epochs.xls for an example of the required format. \n'); end

if ~doit; return; end

%% find files and decide whether to run task;
files=aas_emeg_findfiles(aap,settings.InputFilter,subblock);
if isempty(files); aas_log(aap,1,sprintf('\nFound no data! (Input filter is %s)\n',settings.InputFilter)); end

%% load event specifications from spreadsheet
warning off all; [Numeric,Txt]=xlsread(epochfile); warning on all

%% create sensor*time images (pythagorean sum for grad pairs, so +ive: beware t-test)
S=settings;
S.Fname=char(files);
S.n=147; % appropriate? Odd number better due to sensor positions?
S.pixsize=1; %?
S.interpolate_bad=1;
S.overwrite=settings.Overwrite;
P=spm_eeg_convertmat2ana3D_dm(S);
% option to not overwrite;
% seperately convert mags, grads, eeg, gobal field power for each, and eog;
% small changes for speed and cmdline output
% returns spatial volumes i.e. not eog or gfp
% copes with interpolating channels at edge of array

for f=1:length(files);

    %% load MEG header
    try rehash; load(files{f});
    catch
        try rehash; load(files{f},'-MAT');
        catch; aas_log(aap,1,'\n File may be corrupt? \n');
        end
    end
    fprintf('\nFile: %s',D.fname)
    %% load relevant events/contrasts and time windows of interest
    try
        temptext=regexprep(Txt(1,:),sprintf('%s|',D.events.names{:}),'match');
        col=strcmp(temptext,'match');
    catch continue
    end
    w=Txt(4:5,col);
    TimeWindows=str2num(w{1}); % matrix

    %%  average across each time window for each event
    for tw=1:size(TimeWindows,1);
        twin=TimeWindows(tw,1):TimeWindows(tw,2); %ms
        fprintf('\n\tAveraging across time window %g-%gms',TimeWindows(tw,1),TimeWindows(tw,2))
        swin=floor(twin*D.Radc/1000+D.events.start); % samples
        for et=1:D.Nevents
            Heogvols=strcat(regexprep(files{f},'.mat',''),sprintf('/trialtype%g/average-Heog.img',et));
            Veogvols=regexprep(Heogvols,'Heog','Veog');
            if isempty(findstr(files{f},'-eeg'))
                magvols=regexprep(Heogvols,'Heog','mags');
                gradvols=regexprep(Heogvols,'Heog','grds');
                latvols=regexprep(Heogvols,'Heog','lats');
                longvols=regexprep(Heogvols,'Heog','longs');
                vols=char(magvols,gradvols,latvols,longvols);
                if settings.GFP
                    magvolsRMS=regexprep(Heogvols,'Heog','mags-RMS');
                    gradvolsRMS=regexprep(Heogvols,'Heog','grds-RMS');
                    latvolsRMS=regexprep(Heogvols,'Heog','lats-RMS');
                    longvolsRMS=regexprep(Heogvols,'Heog','longs-RMS');
                    vols=char(vols,magvolsRMS,gradvolsRMS,latvolsRMS,longvolsRMS);
                end
                if settings.EOG
                    vols=char(vols,Heogvols,Veogvols);
                end
            else
                eegvols=regexprep(Heogvols,'Heog','eeg');
                vols=char(eegvols);
                if settings.GFP
                    eegvolsRMS=regexprep(Heogvols,'Heog','eeg-RMS');
                    vols=char(vols,eegvolsRMS);
                end
            end
            vols=cellstr(vols);
            V=spm_vol(vols);
            for v=1:length(V)
                SV=V{v};
                SV.fname=regexprep(V{v}.fname,'average',sprintf('avewin_%g-%gms',min(twin),max(twin)));
                if ~exist(SV.fname,'file') || settings.Overwrite==1
                    try Y=spm_read_vols(V{v});
                    catch
                        if exist(V{v}.fname,'file'); 
                            delete(V{v}.fname)
                            aas_log(aap,true,sprintf('\nFile %s seems corrupt\nIt has been deleted\n',V{v}.fname));
                        else aas_log(aap,true,sprintf('Cannot find file %s\n',V{v}.fname));
                        end
                    end
                    if ~isempty(strfind(V{v}.fname,'eog')) || ~isempty(strfind(V{v}.fname,'RMS'))
                        slice=mean(Y(1,1,swin),3);
                        SV.dim=[1,1,1];
                        SV.mat=[eye(3), ones(3,1); 0 0 0 1];
                    else
                        slice=mean(Y(:,:,swin),3);
                        SV.dim(3)=1;
                        SV.mat(3,4)=1;
                    end
                    try SV=rmfield(SV,{'private','pinfo'}); catch end
                    spm_write_vol(SV,slice);
                end
            end
            fprintf('.')
        end
    end

end % next file

%% smooth sensor-time volumes (from Rik's batch script)
try smooth=str2num(settings.smooth);
catch smooth=[];
end
if ~isempty(smooth)
    for con = 1:size(P,1);
        [pth,nam,ext] = fileparts(P(con,:));
        Pout          = deblank(fullfile(pth,['s' nam ext]));
        if ~exist(Pout,'file') || settings.Overwrite
            fprintf('\nSmoothing %s',P(con,:))
            spm_smooth(spm_vol(P(con,:)),Pout,smooth);
            Pin = strvcat(P(con,:),Pout);
            try spm_imcalc_ui(Pin,Pout,'((i1+eps).*i2)./(i1+eps)',{[],[],'float32',0});   %Reinsert NaNs
            catch
                delete(deblank(P(con,:)));
                error('File %s seems corrupt\nIt has been deleted',deblank(P(con,:)))             
            end
        else fprintf('.')
        end
    end
end

%% mirror, then subtract or average
if settings.mirrorsubtract
    for con = 1:size(P,1);
        [pth,nam,ext] = fileparts(P(con,:));
        if ~isempty(smooth)
            Pin=deblank(fullfile(pth,['s' nam ext]));
            Lout=deblank(fullfile(pth,['ls' nam ext]));
            Bout=deblank(fullfile(pth,['bs' nam ext]));
        else
            Pin=deblank(P(con,:));
            Lout=deblank(fullfile(pth,['l' nam ext]));
            Bout=deblank(fullfile(pth,['b' nam ext]));
        end
        if ~exist(Lout,'file') || settings.Overwrite
            funcstr='i1-flipdim(i1,1)';
            fprintf('\nMirror subtraction of %s',Pin)
            try spm_imcalc_ui(Pin,Lout,funcstr,{[],[],'float32',0});
            catch delete(Pin);
                error('\nFile %s seems corrupt\nIt has been deleted\n',Pin)
            end
        else fprintf('.')
        end
        if ~exist(Bout,'file') || settings.Overwrite
            funcstr='(i1+flipdim(i1,1))/2';
            fprintf('\nBilateral average of %s',Pin)
            try spm_imcalc_ui(Pin,Bout,funcstr,{[],[],'float32',0});
            catch delete(Pin);
                error('\nFile %s seems corrupt\nIt has been deleted\n',Pin)
            end
        else fprintf('.')
        end
    end
end

%             if ~all(cellfun('isempty',regexp(Pout,{'-mags','-longs'})))
%                 funcstr='i1+flipdim(i1,1)';
%                 % symmetric sources produce antisymmetric fields, so sum tests for asymmetry
%                 % (note orientations of longitudinal grads are symmetric)
%             else
%                 funcstr='i1-flipdim(i1,1)';
%                 % symmetric sources produce symmetric fields/potentials, so difference tests for asymmetry
%                 % (note latitudinal grads are asymmetric as they circle the helmet)
%             end

return
