function [aap resp]=aamod_emeg_ploterfs(varargin)
% Plot ERFs/ERPs
% Danny Mitchell 03/03/08

%% check task settings, subject, block etc
[aap subblock doit resp settings]=aa_emeg_checktasksettings(mfilename('fullpath'),varargin);
if ~doit; return; end

%% find files and decide whether to run task;
files=aas_emeg_findfiles(aap,settings.InputFilter,subblock);
if isempty(files); aas_log(aap,1,sprintf('\nFound no data! (Input filter is %s)\n',settings.InputFilter)); end

%% do it

for f=1:length(files);

    % plot timecourses
    [pth xname]=fileparts(files{f});
    if isempty(regexp(files{f},'-eeg','ONCE'))
        if ~exist(fullfile(pth,'figures',[xname '_sensors.ps']),'file') || settings.Overwrite==1;
            aas_emeg_plot_erf(files{f});
        end
    else
        % work out what to do with EEG!
    end

    fprintf('.')
end

return