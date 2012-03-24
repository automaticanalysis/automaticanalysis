function [aap resp]=aamod_emeg_sourcecontrasts(varargin)
% THIS HAD BEEN INCORPORATED INTO AAMOD_EMEG_INVERSION,
% BUT AM NOW USING IT AGAIN SINCE IMPLEMENTING GROUP AND FUSION INVERSIONS
%
% Uses spm_eeg_inv_results to compute induced and/or evoked power of 
% conditional expectation of source activity across time/frequency windows 
% specified in events spreadsheet.
%
% Note: induced contrasts currently do not work for fused data.
%
% (Use aamod_emeg_source2vol to write results to volume, which uses same 
% contrast naming conventions)
%
% Danny Mitchell 02/04/08

%% check task settings, subject, block etc
[aap subblock doit resp settings]=aa_emeg_checktasksettings(mfilename('fullpath'),varargin);

%% CHECK REQUIREMENTS: load spreadsheet containing epoching instructions
epochfile=spm_select('FPList',fullfile(aas_getsesspath(aap,varargin{3},varargin{4}),'events'),'(E|e)pochs.*.xls');
if exist(epochfile,'file')~=2; epochfile=spm_select('FPList',fullfile(aas_getsubjpath(aap,varargin{3}),'events'),'(E|e)pochs.*\.xls'); end
if exist(epochfile,'file')~=2; epochfile=spm_select('FPList',fullfile(aap.acq_details.root,'events'),'(E|e)pochs.*.xls'); end
if exist(epochfile,'file')~=2; error('aa:EpochFileError', '\n Failed to find Epochs.*.xls, required for epoching continuous data. \n This should be placed in "events" folder at either study, subject, or session levels, with lower levels overriding higher ones. \n See examples/Epochs.xls for an example of the required format. \n'); end

if ~doit; return; end

%% find files and decide whether to run task;
files=aas_emeg_findfiles(aap,settings.InputFilter,subblock);
if isempty(files); aas_log(aap,1,sprintf('\nFound no data! (Input filter is %s)\n',settings.InputFilter)); end

%% load event specifications from spreadsheet
warning off all; [Numeric,Txt]=xlsread(epochfile); warning on all

for f=1:length(files);

%% skip EEG or CDA files
    if ~isempty(regexp(files{f},'-eeg','ONCE')); continue; end
    if ~isempty(regexp(files{f},'_CDA_','ONCE')); continue; end
    fprintf('.')
    fprintf('\nFile %s:',files{f});
%% load MEG header
    try rehash; load(files{f});
    catch
        fprintf('\nFailed to load file. Will try again in 10s in case there is an access conflict...\n');
        pause(10)
        try rehash; load(files{f},'-MAT');
        catch fprintf('\nFile may be corrupt.\n'); debugnow
        end
    end

%% load relevant events/contrasts
    try
        % old method:worked when events were labeled from contrast header row
        % temptext=regexprep(Txt(1,:),sprintf('%s|',D.events.names{:}),'match');
        % col=strcmp(temptext,'match');
        
        % new method:works now that events might be labeled from event name column
        anchor=~cellfun('isempty',regexp(D.fname,Txt(1,:)));
        col=~cellfun('isempty',regexp(regexprep(Txt(2,anchor),{'[',']'},{'',''}),regexprep(Txt(2,:),{'[',']'},{'',''})));
    catch continue
    end

%% get windows for time/frequency contrasts
    w=Txt(4:5,col);
    TimeWindows=str2num(w{1}); % matrix
    FrequencyWindows=str2num(w{2}); % matrix

%% find columns containing events and zeroed contrasts
    n=Numeric(:,col);
    n(all(isnan(n')),:)=[]; % remove rows of NaNs; only happened for one subject, possibly due to formatting of .xls
    Ecol=find(size(n,1)-sum(n==0,1)==1); % relevant columns containing events
    Ccol=find(round(sum(n))==0); % relevant columns containing zeroed contrasts
    
%% recode svd values if entered via gui
    if ischar(settings.svd)
        switch settings.svd
            case 'no'; settings.svd=0;
            case 'yes'; settings.svd=1;
            case 'reweighted'; settings.svd=2;
        end
    end
    
%% The penultimate step is to specify a time-frequency window and compute the
        % Conditional expectation of the RMS response.  A simple windowed average
        % is a special case of this, where the frequency is zero.  In this context,
        % the RMS is the same as the absolute value of the time-averaged repsone.  
    tosave=false;
    % for each inversion
    for v=1:length(D.inv);
        D.val=v;
        fprintf('\n - Inversion %g: %s',v,D.inv{v}.comment);
        try D.inv{D.val}.contrast{1}; catch; D.inv{D.val}.contrast={}; end
        % By default only one time/frequency contrast is estimated per inversion,
        % so estimate each contrast on a copy then append to the inversion.
        % Do this both for evoked and induced effects as requested.
        S=D;
        for tw=1:size(TimeWindows,1)
            for fw=1:size(FrequencyWindows,1)
                for induced=0:1
                    if induced && strcmp(settings.type,'evoked'); continue; end
                    if ~induced && strcmp(settings.type,'induced'); continue; end
                    % any other value for settings.type to run both
                    
                    [pth nam]=fileparts(S.fname);
                    S.inv{S.val}.contrast='';
                    warning off all
                    S.inv{S.val}.contrast.woi =TimeWindows(tw,1:2); % peristimulus time (ms)
                    S.inv{S.val}.contrast.fboi=FrequencyWindows(fw,:); % frequency window (Hz)
                    if max(S.inv{S.val}.contrast.fboi)==0; S.inv{S.val}.contrast.fboi=[]; end;
                    warning on all

                    %%%% don't bother breaking up frequencies for short time
                    %%%% window; probably want to (impr|rem)ove this!
                    if ~isempty(S.inv{S.val}.contrast.fboi) && range(S.inv{S.val}.contrast.woi)<100; continue; end
                    %%%%

                    if induced; type='induced'; else type='evoked'; end
                    % could also get evoked effects by running 'induced' on
                    % averaged data
                    
                    if isempty(S.inv{S.val}.contrast.fboi); svdcode=settings.svd;
                    else svdcode=[];
                    end

                    cname=sprintf('%sms_%sHz_%g_%s', ...
                        regexprep(mat2str(S.inv{S.val}.contrast.woi),{' ','[',']'},{'-','',''}), ...
                        regexprep(mat2str(S.inv{S.val}.contrast.fboi),{' ','[',']'},{'-','',''}), ...
                        svdcode,type(1:3));

                    % find contrast if it exists else add new one
                    done=0;
                    for dc=1:length(D.inv{D.val}.contrast)
                        try D.inv{D.val}.contrast{dc};
                            if ~isfield(D.inv{D.val}.contrast{dc},'cname'); break;
                            elseif strcmp(cname,D.inv{D.val}.contrast{dc}.cname); done=1; break;
                            end
                        catch; debugnow % should not happen?
                        end
                        if dc==length(D.inv{D.val}.contrast); dc=dc+1; end
                    end
                    if ~settings.Overwrite && done,
                        fprintf('\nFound contrast window %g: %s, with %g named contrasts', ...
                            dc,cname,length(D.inv{D.val}.contrast{dc}.names));                     
                        continue;
                    else
                        fprintf('\nProcessing %s:',cname);                     
                    end
                    S.inv{S.val}.contrast.type=type;
                    S.inv{S.val}.contrast.svd=settings.svd;
                    S = spm_eeg_inv_results_dm(S); % why does this change S.val?? (e.g. from 1 to 2)
                    % with spm_eeg_inv_results, if using a frequency band,
                    % the fourier set is truncated by svd.
                    % spm_eeg_inv_results_dm gives more options (not sure which
                    % is best!)...
                    % S.svd=0: turn it off
                    % S.svd=1: leave it as it is (to use eigenvectors)
                    % S.svd=2: do it, but then reweight eigenvectors by eigenvalues
                    % (also offers unrecommended rectangular time window, 
                    % and loads data if necessary)

%                     % copy to new field which I'll modify...
%                     S.inv{S.val}.contrast.MW=S.inv{S.val}.contrast.GW;
%                     
%                     % add unsigned (zeroed) contrasts of signed
%                     % localisations (zucosl)
%                     for c=1:length(Ccol)
%                         con=n(:,Ccol(c))'*n(:,Ecol);
%                         uCoSL=zeros(size(S.inv{S.val}.contrast.JW{1}));
%                         for e=1:length(Ecol)
%                             uCoSL=uCoSL+con(e).*S.inv{S.val}.contrast.JW{e};
%                         end
%                         S.inv{S.val}.contrast.MW{end+1}=abs(uCoSL);
%                     end  
%                     
%                     % add signed (zeroed) contrasts of unsigned
%                     % localisations (scoul)
%                     for c=1:length(Ccol)
%                         con=n(:,Ccol(c))'*n(:,Ecol);
%                         sCouL=zeros(size(S.inv{S.val}.contrast.GW{1}));
%                         for e=1:length(Ecol)
%                             sCouL=sCouL+con(e).*S.inv{S.val}.contrast.GW{e};
%                         end
%                         S.inv{S.val}.contrast.MW{end+1}=sCouL;
%                     end 
% 
%                     % zscore power of all contrasts to zero mean and 
%                     % unity standard deviation across mesh
%                     for c=1:length(S.inv{S.val}.contrast.MW)                  
%                         S.inv{S.val}.contrast.MW{c}=zscore(S.inv{S.val}.contrast.MW{c});
%                     end
                    
                    % add contrast names
                    S.inv{S.val}.contrast.names=Txt(1,col);
%                     S.inv{S.val}.contrast.names=[ ...
%                         regexprep(S.inv{S.val}.contrast.names,'.*','$0(zul)') ...
%                         regexprep(S.inv{S.val}.contrast.names(Ccol),'.*','$0(zucosl)') ...
%                         regexprep(S.inv{S.val}.contrast.names(Ccol),'.*','$0(scoul)')];
                    S.inv{S.val}.contrast.cname=cname;

                    % and append contrast to inversion
                    try D.inv{D.val}.contrast{dc}=S.inv{S.val}.contrast;
                    catch D.inv{D.val}.contrast{1}=S.inv{S.val}.contrast;
                    end
                    
                    tosave=true;

                end % evoked then induced
            end % next frequency window
        end % next time window                      
    end % next inversion
    
    if tosave
        fprintf('\nSaving...\n')
        save(fullfile(D.path,D.fname),'D'); % is this best time to save? Can
        % be slow, so don't want to do it too often.
    end
    
end % next file

return
