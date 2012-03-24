function [aap resp]=aamod_emeg_plottopographies(varargin)
% Plot topographies with timepoints in columns, events/contrasts in rows,
% and each data file on a seperate page. 

% Danny Mitchell 03/03/08

%% check task settings, subject, block etc
[aap subblock doit resp settings]=aa_emeg_checktasksettings(mfilename('fullpath'),varargin);
if ~doit; return; end

%% find files and decide whether to run task;
files=aas_emeg_findfiles(aap,settings.InputFilter,subblock);
if isempty(files); aas_log(aap,1,sprintf('\nFound no data! (Input filter is %s)\n',settings.InputFilter)); end

%% do it
    
for f=1:length(files);
    
    load(files{f});
    [pth xname]=fileparts(files{f});
    
    fullwoi=[-D.events.start D.events.stop]*1000/D.Radc; 
    % earliest time, t0, & 6 times after event, with steps multiple of 25ms
    step=ceil((D.events.stop)/7*1000/D.Radc/25)*25;
    snaps=[fullwoi(1) 0:step:fullwoi(2)];
    
    % use the final 0 to set graphics to raster rather than vector
    if ~isempty(regexp(files{f},'-eeg','ONCE'))
        if ~exist(fullfile(pth,'figures',sprintf('%s_topo_EEG.png',xname)),'file') || settings.Overwrite==1;
            aas_emeg_plot_topographies(files{f},'off',snaps,'EEG','A4', ...
                'evts x time',fullfile(pth,'figures',sprintf('%s_topo.png',xname)),0);
        end
    else
        if ~exist(fullfile(pth,'figures',sprintf('%s_topo_Mags.png',xname)),'file') || settings.Overwrite==1;
            aas_emeg_plot_topographies(files{f},'off',snaps,'Mags','A4', ...
                'evts x time',fullfile(pth,'figures',sprintf('%s_topo.png',xname)),0);
        end
        if ~exist(fullfile(pth,'figures',sprintf('%s_topo_Grads.png',xname)),'file') || settings.Overwrite==1;
            aas_emeg_plot_topographies(files{f},'off',snaps,'Grads','A4', ...
                'evts x time',fullfile(pth,'figures',sprintf('%s_topo.png',xname)),0);
        end
    end

end

return