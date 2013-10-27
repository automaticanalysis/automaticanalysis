function EEG = fixegilocs(EEG, fileloc)

alllocs = readlocs(fileloc);

if (EEG.nbchan == 128 || EEG.nbchan == 256) && length(alllocs) == EEG.nbchan+4
    nodatchans = [1 2 3 length(alllocs)]; %3 fiducial and 1 reference channel
    
    alllocs(1).type = 'FID';
    alllocs(2).type = 'FID';
    alllocs(3).type = 'FID';
    alllocs(end).type = 'REF';
    
    EEG.chanlocs = alllocs(setdiff(1:length(alllocs),nodatchans));
    EEG.chaninfo.nodatchans = alllocs(nodatchans);
    EEG.chaninfo.filename = fileloc;
elseif EEG.nbchan == 7 && length(alllocs) == 7
    EEG.chanlocs = alllocs;
else
    error('Mismatch in channel counts: %d in data, %d in %s.', EEG.nbchan, length(alllocs), fileloc);
end