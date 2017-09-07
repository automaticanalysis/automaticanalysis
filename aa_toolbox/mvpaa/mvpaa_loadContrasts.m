% MVPA_LOADCONTRASTS load contrasts for current subject

function contrasts = mvpaa_loadContrasts(aap,p)

% Name of subject...
subjname = aap.acq_details.subjects(p).subjname;

% Get model data from aap
subjmatches=strcmp(subjname,{aap.tasklist.currenttask.settings.model.subject});

% If no exact spec found, try subject wildcard

if (~any(subjmatches))    
    subjwild=strcmp('*',{aap.tasklist.currenttask.settings.model.subject});
    if any(subjwild)
        subjmatches = subjwild;
    end
end

%% Should now have just one model spec
modelnum=find(subjmatches);
if (length(modelnum)>1)
    aas_log(aap,true,sprintf('Error while getting MVPaa contrast details as more than one specification for subject %s',subjname));
end
if (isempty(modelnum))
    aas_log(aap,true,'Cannot find MVPaa contrasts specification. Check either user script');
end

contrasts = aap.tasklist.currenttask.settings.model(modelnum).contrast;