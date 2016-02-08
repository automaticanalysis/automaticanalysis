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
% interest = is this covariate of interest, or a nuisa`  B nce covariate?

function aap=aas_addcovariate(aap,modulename,subject,session,covarName,covarVector,HRF, interest)

if nargin < 7
    HRF = 1; % By default we convolve the covariate
end
if nargin < 8
    interest = 1; % By default, the covariate is of interest
end

% % Get number from end of module name if present in format _%05d (e.g, _00001)
% if (length(modulename>6))
%     moduleindex=str2num(modulename(end-4:end));
%     if (~strcmp(['_' sprintf('%05d',moduleindex)],modulename(length(modulename)-5:end)))
%         moduleindex=1;
%     else
%         modulename=modulename(1:length(modulename)-6);
%     end
% else
%     moduleindex=1;
% end

% Simple way of checking that the covariate is a column vector
if size(covarVector,1) < size(covarVector,2)
    covarVector = covarVector';
end

% Regexp for number at the end of a module name, if present in format _%05d (e.g, _00001)
m1 = regexp(modulename, '_\d{5,5}$');

% Or, we could use '_*' at the end of the module name to specify all modules with that name
m2 = regexp(modulename, '_\*$');

% Or, we might specify certain modules with  '_X/X/X'
m3 = regexp(modulename, '[_/](\d+)', 'tokens');

if ~isempty(m1)
    moduleindex = str2num(modulename(m1+1:end));
    modulename = modulename(1:m1-1);
    
elseif ~isempty(m2)
    modulename = modulename(1:m2-1);
    moduleindex = 1:length(find(strcmp({aap.tasklist.main.module.name}, modulename)));
    
elseif ~isempty(m3)
    modulename = modulename(1:find(modulename=='_',1,'last')-1);
    moduleindex = cellfun(@str2num, [m3{:}]);
    
else
    moduleindex = 1;
end

% find model that corresponds and add contrast to this if it exists
for m = 1 : length(moduleindex)
    
    mInd = moduleindex(m);
    
    % find model that corresponds and add event to this if it exists
    whichmodel=[strcmp({aap.tasksettings.(modulename)(mInd).modelC.subject},subject)] & ...
        [strcmp({aap.tasksettings.(modulename)(mInd).modelC.session},session)];
    
    if (~any(whichmodel))
        emptymod=[];
        emptymod.subject=subject;
        emptymod.session=session;
        emptymod.covariate.name=covarName;
        emptymod.covariate.vector=covarVector;
        emptymod.covariate.HRF=HRF;
        emptymod.covariate.interest=interest;
        aap.tasksettings.(modulename)(mInd).modelC(end+1)=emptymod;
    else
        aap.tasksettings.(modulename)(mInd).modelC(whichmodel).covariate(end+1).name=covarName;
        aap.tasksettings.(modulename)(mInd).modelC(whichmodel).covariate(end).vector=covarVector;
        aap.tasksettings.(modulename)(mInd).modelC(whichmodel).covariate(end).HRF=HRF;
        aap.tasksettings.(modulename)(mInd).modelC(whichmodel).covariate(end).interest=interest;
    end
end