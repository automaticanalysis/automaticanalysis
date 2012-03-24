function [aap resp]=aamod_emeg_epoch(varargin)
% Epochs EMEG data from spreadsheet containing event definitions. This
% should be called Epochs*.xls and can be placed in events folder at study,
% subject or session level; lower levels override higher ones. For an
% example of the required format, see examples/Epochs.xls.
%
% (Older version also had settings for resychronising event timings and
% only epoching selected events, but now specify these in Epochs*.xls.)
%
% Danny Mitchell 03/03/08

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
Numeric(all(isnan(Numeric')),:)=[];

%% run task for each file
for f=1:length(files);
    
%% prepare parameters
    E.events.Inewlist=0; % probably no need to recode triggers?
    E.D=files{f};
    [pth name]=fileparts(E.D); 
    name=strrep(name,'_raw',''); % remove raw suffix from continuous data
       
%% loop through columns skipping those where a subsequent column contains the same events
    for c=4:size(Numeric,2)
        % check for redundancy...
        eventindices=find(Numeric(:,c)~=0); % events in this column
        temp=Numeric(:,4:end); % 1st three columns have trigger codes, event names and resynch options
        if sum(all(temp(eventindices,:)~=0,1))>1
            Numeric(:,c)=zeros(size(Numeric,1),1);
            continue;
        end

%% epoch all nonzero events in column, using window size and resynch specified in Epochs.xls.
        outfile=fullfile(pth,['e_' name '_' Txt{1,c} '.mat']);
        if ~exist(outfile,'file') || settings.Overwrite==1;
            E.events.types=Numeric(eventindices,1);
            if isempty(E.events.types); continue; end
            epochwindow=str2num(Txt{2,c}); % vector
            E.events.start=epochwindow(1);
            E.events.stop=epochwindow(2);
            E.events.resynch=Numeric(eventindices,2);
            fprintf('\nEpoching %s from %s, using %s...',Txt{1,c},name,epochfile)
            D=spm_eeg_epochs(E); % version in cbu_updates can cope with resynching different event types
            
%% warn and delete file if it has no events (hopefully shouldn't happen!), 
%% otherwise add any ROI information and save with column heading as suffix
            if D.events.Ntypes>0
                % Rename according to final contrast
                newfname=['e_' name '_' Txt{1,c} '.mat'];
                newfnamedat=['e_' name '_' Txt{1,c} '.dat'];
                save(fullfile(D.path, D.fname), 'D');
                cmd=['mv ' fullfile(D.path,D.fname) ' ' fullfile(D.path,newfname)]; unix(cmd);
                cmd=['mv ' fullfile(D.path,D.fnamedat) ' ' fullfile(D.path,newfnamedat)]; unix(cmd);
                D.fname=newfname;
                D.fnamedat=newfnamedat;
                % add event names
                if length(eventindices)~=length(D.events.types); error('\nCheck epoch file.\n'); end
                D.events.names=Txt(7+eventindices,3)';
                % add ROI to structure if necessary
                try
                    D.ROIs.space=Txt(6,c);
                    D.ROIs.co=str2num(char(Txt(7,c))); % vector
                catch
                end
                save(fullfile(D.path, D.fname), 'D');

            else
                fprintf('\n - Warning: Found no events for %s: Deleting.',Txt{1,c});
                delete(fullfile(D.path,D.fnamedat));
                delete(fullfile(D.path,D.fname));
                continue
            end
        end
    end % next event/contrast column 

    fprintf('.')
end % next fif file

return

%% some old code for epoching from MaxAve parameter file 
%         % epoch all event types specified using parameters saved from MaxAve
%         % THIS IS OLD - NOT SURE IF IT STILL WORKS.
%         if exist(fullfile(aapMEG.AnalysisDirectory,aap.acq_details.subjects(subblock{1}).Name,aap.acq_details.sessions(subblock{2}).Name,['e_' nam '.dat']),'file')==0 || settings.Overwrite==1;
%             fprintf('\nEpoching from MaxAve parameter file...')
%             try
%                 E.events.types=D.events.types;
%             catch
%                 load(strrep(outfile,'.dat','.mat'));
%                 E.events.types=D.events.types;
%             end
%             E.twin=settings.twin;
%             % extend time window to include time windows of all specified
%             % events
%             if E.twin(1)==-Inf
%                 load(strrep(infile,'.fif','_ave_nave'),'windows'); % load time windows or each average from mat file
%                 E.twin(1)=min(windows(:,1));
%             end
%             if E.twin(2)==Inf
%                 load(strrep(infile,'.fif','_ave_nave'),'windows'); % load time windows or each average from mat file
%                 E.twin(2)=max(windows(:,2));
%             end
%             E.events.start=E.twin(1);
%             E.events.stop=E.twin(2);
%             spm_eeg_epochs(E);
%         end