% Allows for specifying contrast vector in an "easy" way:
% One needs to specify contrast weight only for Event-of-Interest. Contrast vector will be filled with zeroes (e.g. for Event for erroneous scans).
% Tibor Auer MRC CBU Cambridge 2012-2013

function con = aas_contrast2model(aap,nModel,convect)
nEV = numel(aap.tasksettings.aamod_firstlevel_model.model(nModel).event);
con = zeros(1,nEV);
con(1:numel(convect)) = convect;