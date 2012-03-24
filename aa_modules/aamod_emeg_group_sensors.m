function [aap resp]=aamod_emeg_group_sensors(varargin)
% Danny Mitchell 13/03/08; parallelised 11/08

%% check task settings, subject, block etc
[aap subblock doit resp settings]=aa_emeg_checktasksettings(mfilename('fullpath'),varargin);  

%% CHECK REQUIREMENTS: load spreadsheet containing epoching instructions
epochfile=spm_select('FPList',fullfile(aap.acq_details.root,'events'),'(E|e)pochs.*.xls');
if exist(epochfile,'file')~=2; error('aa:EpochFileError', '\n Failed to find Epochs.*.xls, required for epoching continuous data. \n This should be placed in "events" folder at either study, subject, or session levels, with lower levels overriding higher ones. \n See examples/Epochs.xls for an example of the required format. \n'); end

if ~doit; return; end

%%
gfile=fullfile(aap.acq_details.root,...
    sprintf('%s_%g.parallel.mat',mfilename,(aap.tasklist.currenttask.index)));

%%
if (isempty(varargin) && ~exist(gfile,'file')) || (nargin>1 && strcmp(varargin{2},'parallelise'))
    % prepare globals and looping variable
    %% get files to average
    fprintf('\nFinding files from %g blocks of %g subjects', ...
        length(aap.acq_details.sessions),length(aap.acq_details.subjects))

    erfs=[];eblocks=[];
    settings.InputFilter=strrep(settings.InputFilter,'BLOCK','');

    for s=1:length(aap.acq_details.subjects)
        fprintf('.')
        for b=aap.acq_details.selected_sessions
            block=aap.acq_details.sessions(b).name;
            subdir=fullfile(aap.acq_details.root,aap.acq_details.subjects(s).megname,block);

            % get all erf files
            temp=spm_select('List',subdir,settings.InputFilter); % not tf data
            erfs=char(erfs,temp);
            eblocks=char(eblocks,repmat(block,[size(temp,1) 1]));

        end % next block
    end % next subject

    try erfs(1,:)=[]; catch end
    try eblocks(1,:)=[]; catch end

    % get unique names
    if isempty(erfs);
        aas_log(aap,true,sprintf('Found no data! Input filter is: %s', ...
            settings.InputFilter));
    else
        [enames I]=unique(erfs,'rows');
        eblocks=eblocks(I,:);
    end

    %% fill (contrast x subject) matrix of full file paths
    efiles={};
    warnings={};
    for s=1:length(aap.acq_details.subjects)
        for e=1:size(enames,1)
            filename=deblank(fullfile(aap.acq_details.root,aap.acq_details.subjects(s).megname,eblocks(e,:),enames(e,:)));
            if exist(filename,'file');
                efiles(e,s)=cellstr(filename);

                %             %%% debugging: check that no files are corrupt
                %             try load(filename)
                %             catch warnings=[warnings;filename];
                %             end
                %             %%%

            else efiles(e,s)={' '};
            end
        end
    end
    clear eblocks erfs

    %%% debugging: check that no files are corrupt
    fprintf('\nSummary of subject numbers per event:\n')
    [enames num2str(sum(cellfun('isempty',regexp(efiles,'^ $')),2))]
    if ~isempty(warnings)
        fprintf('\nPossibly corrupt files! Please check:\n')
        char(warnings)
        keyboard
    end
    %%%

    %% confirm/create group analysis directory
    try cd(fullfile(aap.acq_details.root,'GroupAnalysis_Sensors'));
    catch
        mkdir(fullfile(aap.acq_details.root,'GroupAnalysis_Sensors'));
        cd(fullfile(aap.acq_details.root,'GroupAnalysis_Sensors'));
    end
    if ~exist('figures','dir'); mkdir('figures'); end

    %% save variables for parallelisation
    loopvar=cellstr(enames);
    save(gfile,'efiles','epochfile','loopvar');
    %aap.tasksettings.(mfilename).PARALLEL=gfile;

    if ~isempty(varargin); return; end
end

%% do it
try load(gfile); loopvar;
catch aas_log(aap,1,sprintf('\n***Failed to prepare parallelisation structure!\n%s',resp));
end;

if isempty(subblock); subblock=1:length(loopvar); 
else subblock=subblock{1};
end
enames=char(loopvar);

for e=subblock
    
    %% load event specifications from spreadsheet
    warning off all; [Numeric,Txt]=xlsread(epochfile); warning on all

    ename=deblank(enames(e,:));
    [pth nam ext]=fileparts(ename);

    S.P=char(efiles{e,:}); % contrast files from each subject
    S.P(all(S.P==' ',2),:)=[]; % remove missing files
    if size(S.P,1)<2;  % only continue if there are multiple subjects
        fprintf('\nFailed to find multiple subjects for contrast %s.', ename)
        continue;
    else fprintf('\nProcessing %s (for %g subjects)...', ename, size(S.P,1))
    end

    % load(deblank(S.P(1,:))); % load data structure from any subject

    %%% calculate grandmean
    Gmean=fullfile(aap.acq_details.root,'GroupAnalysis_Sensors',[sprintf('g%g',size(S.P,1)) nam ext]);
    if ~exist(Gmean,'file') || settings.Overwrite==1
        % delete any old averages with fewer subjects
        %             old=cellstr(spm_select('List',pwd,nam));
        %             for o=1:length(old)
        %                 warning off all; delete(old{o}); warning on all
        %             end
        S.Pout=Gmean;
        D=spm_eeg_grandmean(S); % average across subjects
        try D=rmfield(D,{'inv','val'}); catch; end;  % remove any inversions
        save(fullfile(D.path,D.fname),'D');
    else
        load(Gmean);
    end
    
    if isfield(D,'Nfrequencies'); continue; end

%% load time-windows of interest
    temptext=regexprep(Txt(1,:),sprintf('%s|',D.events.names{:}),'match');
    fc=find(strcmp(temptext,'match'));
    woi=str2num(Txt{4,fc(1)}); % matrix
    woi=woi(end,:); % just use later time window
    
%% create group summaries of field topographies
    % load data structure from any subject
    %if ~exist('D','var'); load(deblank(S.P(1,:))); end

    % (average of woi, subs x subs, in page layout)
    if ~isempty(regexp(D.fname,'-eeg','ONCE'))
        outfile=fullfile(aap.acq_details.root,'GroupAnalysis_Sensors','figures',sprintf('g%g%s_topo_EEG_%02g.png',size(S.P,1),nam, D.Nevents));
        if ~exist(outfile,'file') || settings.Overwrite==1;
            aas_emeg_plot_topographies(char(S.P,Gmean),'off',woi,'EEG','A4','subs x subs',outfile,0);
        end
    else
        outfile=fullfile(aap.acq_details.root,'GroupAnalysis_Sensors','figures',sprintf('g%g%s_topo_Mags_%02g.png',size(S.P,1),nam, D.Nevents));
        if ~exist(outfile,'file') || settings.Overwrite==1;
            aas_emeg_plot_topographies(char(S.P,Gmean),'off',woi,'Mags','A4','subs x subs',outfile,0);
        end
        outfile=fullfile(aap.acq_details.root,'GroupAnalysis_Sensors','figures',sprintf('g%g%s_topo_Grads_%02g.png',size(S.P,1),nam, D.Nevents));
        if ~exist(outfile,'file') || settings.Overwrite==1;
            aas_emeg_plot_topographies(char(S.P,Gmean),'off',woi,'Grads','A4','subs x subs',outfile,0);
        end
    end

%% video of grandmean, evts x evts, in screen layout, defining colormap
% from mean of 50-200ms
%     if D.Nsamples<300; msperframe=4; else msperframe=8; end % might need to conserve memory
%     if ~isempty(regexp(D.fname,'-eeg','ONCE'))
%         outfile=fullfile(aap.acq_details.root,'GroupAnalysis_Sensors','figures',[sprintf('g%g',size(S.P,1)) nam '_topo_EEG.avi']);
%         if ~exist(outfile,'file') || settings.Overwrite==1
%             aas_emeg_plot_topographies(S.Pout,msperframe,[50 200],'EEG','Screen','evts x evts',outfile,0);
%         end
%     else
%         outfile=fullfile(aap.acq_details.root,'GroupAnalysis_Sensors','figures',[sprintf('g%g',size(S.P,1)) nam '_topo_Mags.avi']);
%         if ~exist(outfile,'file') || settings.Overwrite==1
%             aas_emeg_plot_topographies(S.Pout,msperframe,[50 200],'Mags','Screen','evts x evts',outfile,0);
%         end
%         outfile=fullfile(aap.acq_details.root,'GroupAnalysis_Sensors','figures',[sprintf('g%g',size(S.P,1)) nam '_topo_Grads.avi']);
%         if ~exist(outfile,'file') || settings.Overwrite==1
%             aas_emeg_plot_topographies(S.Pout,msperframe,[50 200],'Grads','Screen','evts x evts',outfile,0);
%         end
%     end

%% plot grandmean ERFs?
    [path name]=fileparts(Gmean);
    outfile=fullfile(path,'figures',[name '_sensors.ps']);
    if ~exist(outfile,'file') || settings.Overwrite==1;
        try 
            aas_emeg_plot_erf(Gmean);
        catch
            % keyboard % work out what to do for EEG!
            fprintf('\nNot plotting ERPs for %s',S.Pout');
        end
    end
    
    %%% create concatenation of subject ERFs for ease of plotting group?
    %% would be better not to accumulate bad channels, but to interpolate
    %% them before merging?
    S.Pout=fullfile(aap.acq_details.root,'GroupAnalysis_Sensors',[sprintf('c%g',size(S.P,1)) nam ext]);
    if ~exist(S.Pout,'file') || settings.Overwrite==1
        S.D=S.P;
        % don't want to recode, but need to specify this anyway
        if ~exist('D','var'); load(deblank(S.P(1,:))); end
        S.recode=repmat({D.events.types},[size(S.P,1) 1]);
        D=spm_eeg_merge_dm(S); % allow specification of output file
        try D=rmfield(D,{'inv','val'}); catch end; % remove any inversions
        save(fullfile(D.path,D.fname),'D');
    end
    
    fprintf('.')
end
 
return