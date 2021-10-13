% Adds a model for group stat
% FORMAT function aap = aas_add_meeg_groupmodel(aap, modulename, subjspec, trialmodelspec, channelspec, groupmodelspec, timewindowspec, modelname)
%   - subjspec: cell array of subjectnames to be included in the model
%   - trialmodelspec: cell array of trialmodels or peaks (see aamod_meeg_timelockanalysis) for each participant is included in the dependent variable. For
%       peaks, the cell array should contain a trialmodelspec and a peakspec delimited with "_" (e.g. "STIM_P300")
%       N.B.: The vector of dependent variable is compiled by adding all trialmodels for participant #1, then the same for all participants.
%   - channelspec: describing how channels should be considered. It can be
%       - 'all': when channels provide topographic information
%       - cell array of <channel label>: when data only from the selected channels are considered
%       - cell array of <Nxchannel label>s: when each channel is combined based on its weight N, and the combined activity is included in the dependent variable. 
%   - groupmodelspec: model describing the elements of the independent variable
%       N.B.: It must be a vector of 1xnumber of dependent element, i.e. subjects*trialmodels*channels(if specified as cell array)
%       N.B.: When entering two samples, it specifies a T-test with groupmodelspec==1 minus groupmodelspec==2
%   - timewindowspec: timewindow of interest specified as [start stop] in millisecond or 'all'
%   - modelname: name of the model (used for display and saving)
% N.B.: '*' for subjspec or trialmodelspec means all subjects or trialmodels, respectively
%
% Examples
%aap = aas_add_meeg_groupmodel(aap, 'aamod_meeg_timelockstatistics', '*', '*', repmat([1,1],1,12), 'all');
%aap = aas_add_meeg_groupmodel(aap, 'aamod_meeg_timelockstatistics', {'003' '004' '005'}, 'STIM', ones(1,3), 'all');

function aap = aas_add_meeg_groupmodel(aap, modulename, subjspec, trialmodelspec, channelspec, groupmodelspec, timewindowspec, modelname)

% Regexp for number at the end of a module name, if present in format _%05d (e.g. _00001)
m1 = regexp(modulename, '_\d{5,5}$');

% Or, we could use '_*' at the end of the module name to specify all modules with that name
m2 = regexp(modulename, '_\*$');

% Or, we might specify certain modules with  '_X/X/X' (e.g. _00001/00002/00004)
m3 = regexp(modulename, '[_/](\d+)', 'tokens');

if ~isempty(m1)
    moduleindex = str2double(modulename(m1+1:end));
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

if ~exist('modelname','var')
    modelname='';
end
if any(modelname == '_'), aas_log(aap,true,sprintf('ERROR: Modelname %s has "_"!',modelname)); end

if ~iscell(subjspec), subjspec = {subjspec}; end
if subjspec{1} == '*'
    subjspec = {aap.acq_details.subjects.subjname};
end

if iscell(channelspec) && (any(channelspec{1}== '+') || any(channelspec{1}== '-'))
    weights = cellfun(@(x) str2double(x), regexp(channelspec,'[\+\-0-9]*','match'));
    channels = regexp(channelspec,'(?<=x)[a-zA-Z0-9]*','match','once');
    channelspec = [];
    channelspec.channels = channels;
    channelspec.weights = weights;
end

% find model that corresponds and add trialmodel to this if it exists
for m = 1:numel(moduleindex)
    
    mInd = moduleindex(m);
    
    if any(strcmp({aap.tasksettings.(modulename)(mInd).model.name},modelname))
        aas_log(aap,true,['Model ' modelname ' has been already specified.']);
    end
    
    emptymod = aap.tasksettings.(modulename)(mInd).model(1); % The first one is usually empty, makes for a good template in case the structure changes
    emptymod.name = modelname;
    emptymod.subjects = subjspec;
    if ~iscell(trialmodelspec), trialmodelspec = {trialmodelspec}; end
    if trialmodelspec{1} == '*'
        % specify trialmodels based on the first subject (in the model)
        trialmodelspec = strcmp({aap.tasksettings.aamod_meeg_timelockanalysis(mInd).trialmodel.subject},emptymod.subjects{1});
        trialmodelspec = [aap.tasksettings.aamod_meeg_timelockanalysis(mInd).trialmodel(trialmodelspec).model];
        trialmodelspec = {trialmodelspec.name};
    end
    emptymod.trialmodel = trialmodelspec;
    if size(groupmodelspec,1) > size(groupmodelspec,2)
        aas_log(aap,false,sprintf('Groupmodel specification MUST be a 1xN matrix. The entered matrix with size %dx%d will be transposed',size(groupmodelspec)));
        groupmodelspec = groupmodelspec';
    end
    if numel(groupmodelspec) ~= (numel(emptymod.subjects) * numel(emptymod.trialmodel))
        aas_log(aap,true,sprintf('Number of elements in groupmodel (%d) does not match with number of subjects (%d) x number of trialmodels (%d).',...
            numel(groupmodelspec),numel(emptymod.subjects),numel(emptymod.trialmodel)));
    end
    emptymod.channels = channelspec;
    emptymod.timewindow = timewindowspec;
    emptymod.groupmodel = groupmodelspec;

    aap.tasksettings.(modulename)(mInd).model(end+1)=emptymod;
end
