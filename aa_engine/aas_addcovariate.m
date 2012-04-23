% Adds a covariate to a model
% function
% aap=aas_addcovariate(aap,modulename,subject,session,covarName,covariate, HRF, interest)
% 
% modulename=name of module (e.g.,'aamod_firstlevel_model') for which this
%   covariate applies
% subject = subject for whom this model applies
% session = session for which this applies
% covarName = name of the covariate
% covarVector = covariate vector, which should be as long as the session
% HRF = do we want to convolve this covariate with the HRF? (0 - no; 1 - yes)
% interest = is this covariate of interest, or a nuisance covariate?

function aap=aas_addcovariate(aap,modulename,subject,session,covarName,covarVector,HRF, interest)

if nargin < 7
    HRF = 1; % By default we convolve the covariate
end
if nargin < 8
    interest = 1; % By default, the covariate is of interest
end

% Get number from end of module name if present in format _%05d (e.g, _00001)
if (length(modulename>6))
    moduleindex=str2num(modulename(end-4:end));
    if (~strcmp(['_' sprintf('%05d',moduleindex)],modulename(length(modulename)-5:end)))
        moduleindex=1;
    else
        modulename=modulename(1:length(modulename)-6);
    end
else
    moduleindex=1;
end

% find model that corresponds and add event to this if it exists
whichmodel=[strcmp({aap.tasksettings.(modulename)(moduleindex).modelC.subject},subject)] & ...
    [strcmp({aap.tasksettings.(modulename)(moduleindex).modelC.session},session)];

if (~any(whichmodel))
    emptymod=[];
    emptymod.subject=subject;
    emptymod.session=session;
    emptymod.covariate.name=covarName;
    emptymod.covariate.vector=covarVector;
    emptymod.covariate.HRF=HRF;
    emptymod.covariate.interest=interest;
    aap.tasksettings.(modulename)(moduleindex).modelC(end+1)=emptymod;
else
    aap.tasksettings.(modulename)(moduleindex).modelC(whichmodel).covariate(end+1).name=covarName;
    aap.tasksettings.(modulename)(moduleindex).modelC(whichmodel).covariate(end).vector=covarVector;
    aap.tasksettings.(modulename)(moduleindex).modelC(whichmodel).covariate(end).HRF=HRF;
    aap.tasksettings.(modulename)(moduleindex).modelC(whichmodel).covariate(end).interest=interest;
end