function meeg_diagnostics_continuous(EEG,diag,figtitle,savepath)
if nargin < 3, figtitle = 'Sample'; end
if ~isempty(diag.freqrange)
    f = figure('Name',figtitle);
    if isempty(diag.freq)
        pop_spectopo(EEG,1,[0 650399],'EEG','freqrange',diag.freqrange,'electrodes','off');
    else
        pop_spectopo(EEG,1,[0 650399],'EEG','freq',diag.freq,'freqrange',diag.freqrange,'electrodes','off');
    end
    if nargin >= 4
        set(f,'Position',[1 1 1920 1080]);
        set(f,'PaperPositionMode','auto');
        print(f,'-noui',savepath,'-djpeg','-r150');
        close(f);
    end
end