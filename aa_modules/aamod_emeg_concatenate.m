function [aap resp]=aamod_emeg_concatenate(varargin)
% Epoched data from sequentially numbered blocks of the same type are 
% combined using spm_eeg_merge. This just concatenates the files, 
% preserving all the epochs.
%
% Danny Mitchell 03/03/08

%% check task settings, subject, block etc
[aap subblock doit resp settings]=aa_emeg_checktasksettings(mfilename('fullpath'),varargin);
if ~doit; return; end

%% find files and decide whether to run task;
files=aas_emeg_findfiles(aap,settings.InputFilter,subblock);
if isempty(files); aas_log(aap,1,sprintf('\nFound no data! (Input filter is %s)\n',settings.InputFilter)); end

%% find filter for each group of files to be averaged; assume repeated blocks are numbered
q=regexprep(files,{'(.*/)',sprintf('(%s)(\\d)',aap.acq_details.sessions(subblock{2}).name)},{'^','$1[\\d]*'});
groups={};
for candidate=1:numel(q)
    if sum(strcmp(groups,q{candidate}))==0
        groups=[groups; q{candidate}];
    end
end

pth=fileparts(files{1});
%% combine spm files within each group
for e=1:length(groups)
    S.D=spm_select('List',pth,groups{e});
    fprintf('.');
    if ~exist([pth filesep 'c' deblank(S.D(1,:))],'file') || settings.Overwrite==1;
        if size(S.D,1)>1
            % merge
            S.D=[repmat([pth filesep],size(S.D,1),1) S.D];
            S.recode={};
            for f=1:size(S.D,1)
                % don't want to recode, but need to specify this anyway
                load(deblank(S.D(f,:)));
                S.recode=[S.recode; D.events.types];
            end
            spm_eeg_merge(S);
        else
            % single file; nothing to merge; just appeand c to name
            load([pth filesep deblank(S.D(1,:))]);
            D.fname=['c' D.fname];
            save(fullfile(D.path, D.fname), 'D');
        end
    else continue % next epoch
    end
end

return