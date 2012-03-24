function [aap resp]=aamod_emeg_splitsensors(varargin)
% Split MEG data into seperate files for magnetometers and gradiometers,
% using Rik's spm_eeg_splitFIF.
%
% Danny Mitchell 03/03/08

%% check task settings, subject, block etc
[aap subblock doit resp settings]=aa_emeg_checktasksettings(mfilename('fullpath'),varargin);
if ~doit; return; end

%% find files and decide whether to run task;
files=aas_emeg_findfiles(aap,settings.InputFilter,subblock);
if isempty(files); aas_log(aap,1,sprintf('\nFound no data! (Input filter is %s)\n',settings.InputFilter)); end

%% check whether each file needs splitting or has already been split
for f=1:length(files);  
    if ~isempty(regexp(files{f},'(-mags|-grds|-eeg)','ONCE')); continue; end
    [pth nam ext]=fileparts(files{f});
    if ~settings.Overwrite && ...
            exist(fullfile(pth,['s' nam '-mags' ext]),'file') && ...
            exist(fullfile(pth,['s' nam '-grds' ext]),'file');
        fprintf('.')
        continue;
    end
    
%% split file
    S.D=files{f};
    S.grms=false;
    D=spm_eeg_splitFIF(S);
    
%% add 's' prefix to outputs
    rehash
    movefile(fullfile(D.path,D.fname),fullfile(D.path,['s' D.fname]));
    movefile(fullfile(D.path,D.fnamedat),fullfile(D.path,['s' D.fnamedat]));
    movefile(fullfile(D.path,strrep(D.fname,'grds','mags')),fullfile(D.path,['s' strrep(D.fname,'grds','mags')]));
    movefile(fullfile(D.path,strrep(D.fnamedat,'grds','mags')),fullfile(D.path,['s' strrep(D.fnamedat,'grds','mags')]));
    D.fname=['s' D.fname];
    D.fnamedat=['s' D.fnamedat];
    save(fullfile(D.path, D.fname), 'D');
    load(fullfile(D.path,strrep(D.fname,'grds','mags')));
     D.fname=['s' D.fname];
    D.fnamedat=['s' D.fnamedat];
    save(fullfile(D.path, D.fname), 'D');   
end % next file

return

