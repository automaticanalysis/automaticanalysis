function [d, seriesnum] = aas_get_series(aap,modality,subj,sess)

switch modality
    case {'FMRI' 'functional'}
        seriesnumbers = aap.acq_details.subjects(subj).seriesnumbers;
    case {'MEG' 'EEG' 'meg'}
        seriesnumbers = aap.acq_details.subjects(subj).megseriesnumbers;
    case {'DWI' 'diffusion'}
        seriesnumbers = aap.acq_details.subjects(subj).diffusion_seriesnumbers;
    case {'MTI' 'special'}
        seriesnumbers = aap.acq_details.subjects(subj).specialseries;
end


outsess = sess;
for d = 1:numel(aap.acq_details.subjects(subj).mriname)
    series = seriesnumbers{d};
    outsess = outsess - numel(series);
    if outsess < 1
        outsess = numel(series) + outsess;
        break;
    end
end

if iscell(series) % Multi-meas (e.g. multi-echo)
    seriesnum=series{outsess};
else
    seriesnum=series(outsess);
end