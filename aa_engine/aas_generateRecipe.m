%% aas_generateRecipe(aap, analyses)
function currentRecipe = aas_generateRecipe(rootPth, analyses, prefix)

if nargin < 3
    prefix = 'aap_recipe_';
end

if ~isunix
    error('This function (currently) only runs on Unix based systems currently.')
end

% Need to add path to recipe so that it can be detected...
if ~exist(rootPth, 'dir')
    mkdir(rootPth)
end
addpath(rootPth)

currentRecipe = fullfile(rootPth, 'aap_tasklist_current.xml');

% Remove current tasklist file and create empty one
try
    unix(['rm ' currentRecipe]);
catch
end
unix(['touch ' currentRecipe]);

% Document begin and initialisation modules
unix(['cat ' which('aap_tasklist_begin.xml') ' >> ' currentRecipe]);

% Analyses modules
for iAnalysis = 1:length(analyses)
    moduleAdded = which([prefix analyses{iAnalysis} '.xml']);
    if isempty(moduleAdded)
        error(sprintf('Cannot find tasklist %s ', [prefix analyses{iAnalysis} '.xml']))
    end
    unix(['cat ' moduleAdded ' >> ' currentRecipe]);    
end

% Garbage collection module and document end
unix(['cat ' which('aap_tasklist_end.xml') ' >> ' currentRecipe]);

% Make file executable and group private...
unix(['chmod 770 ' currentRecipe]);