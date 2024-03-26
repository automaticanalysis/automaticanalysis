function [d, seriesnum] = aas_get_series(aap,modality,subj,sess)

switch modality
    case {'FMRI' 'functional' 'session'}
        seriesnumbers = aap.acq_details.subjects(subj).seriesnumbers;
    case {'MEEG' 'MEG' 'EEG' 'meg'}
        seriesnumbers = aap.acq_details.subjects(subj).meegseriesnumbers;
    case {'DWI' 'diffusion'}
        seriesnumbers = aap.acq_details.subjects(subj).diffusion_seriesnumbers;
    case {'MTI' 'special'}
        seriesnumbers = aap.acq_details.subjects(subj).specialseries;
end


outsess = sess; hasSess = false;
for d = 1:numel(aap.acq_details.subjects(subj).mriname)
    series = seriesnumbers{d};
    outsess = outsess - numel(series);
    if outsess < 1
        hasSess = true;
        outsess = numel(series) + outsess;
        break;
    end
end

if hasSess
    if iscell(series) % Multi-meas (e.g. multi-echo)
        seriesnum=series{outsess};
    else
        seriesnum=series(outsess);
    end
else
    d = [];
    seriesnum = [];
end