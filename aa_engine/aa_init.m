% Automatic analysis - initialise parameters from recipe
% Initialise parameters using recipe by loading all of the parameter sets
% that the recipe specifies.
% Rhodri Cusack 2004-6 MRC CBU Cambridge

function [aap]=aap_init(aap)

global defaults

% Set SPM defaults
spm_defaults

aap.spm.defaults=defaults;

% Load defaults as specified in recipe
aap=feval(aap.recipe.directory_conventions,aap);
aap=feval(aap.recipe.options,aap);
aap=feval(aap.recipe.spmanalysis,aap);
aap=feval(aap.recipe.acq_details,aap);
aap=feval(aap.recipe.tasklist,aap);
try
    aap=feval(aap.recipe.meg,aap);
catch
end;

if (isfield(aap.recipe,'special'))
    aap=feval(aap.recipe.special,aap);
end;

% store a copy of the virgin recipe before the user screws it up
aap.aap_beforeuserchanges=[];
aap.aap_beforeuserchanges=aap;
