function [aap resp]=aamod_emeg_filter(varargin)
% Apply temporal filters using spm_eeg_filter. This uses a 5th order
% butterworth filter.
% Shouldn't it be possible to apply mulitple filters at once instead of
% sequentially?
% Use modified version of spm_eeg_filter (was spm_eeg_filter_dm) in
% cbu_updates as it's faster
% Danny Mitchell 04/04/08

%% check task settings, subject, block etc
[aap subblock doit resp settings]=aa_emeg_checktasksettings(mfilename('fullpath'),varargin);
if ~doit; return; end

%% find files and decide whether to run task;
files=aas_emeg_findfiles(aap,settings.InputFilter,subblock);
if isempty(files); aas_log(aap,1,sprintf('\nFound no data! (Input filter is %s)\n',settings.InputFilter)); end

%% run task for each file

for f=1:length(files);
    [pth name ext]=fileparts(files{f});
    outputfile=fullfile(pth,['f' name ext]);
    fprintf('.')
        
    % output file is created before filtering, but D.fname is only updated
    % afterwards, so check this to see if completed
    try load(outputfile)
        if strcmp(D.fname(1),'f') && ~settings.Overwrite; continue; end
    catch % do it
    end
    
    S.filter=settings;
    S.D = files{f};
    S.filter.band = settings.band;
    S.filter.PHz = str2num(settings.PHz); % scalor or 2 vector
    tic;
    spm_eeg_filter(S); % modified so less file access, so hopefully much quicker; now in cbu_updates
    %fprintf('Took %g minutes',toc/60);

    % no longer used...
%     if filt<length(settings.band) % if want to apply another filter...
%         % remove 'f' prefix ready for next filter
%         load(outputfile)
%         cmd=['mv ' fullfile(pth,['f' name '.mat']) ' ' fullfile(pth,[name '.mat'])]; unix(cmd);
%         cmd=['mv ' fullfile(pth,['f' name '.dat']) ' ' fullfile(pth,[name '.dat'])]; unix(cmd);
%         D.fname=[name '.mat'];
%         D.fnamedat=[name '.dat'];
%         save(fullfile(D.path, D.fname), 'D');
%     end  

end

return