% Adds a definition for peak within a trial
% FORMAT function aap = aas_add_meeg_peakdef(aap, modulename, subjname, trialname, direction, toi, peakname)
%   - aap: aap structure with parameters and tasklist
%   - modulename: name of module (e.g.,'aamod_meeg_timelockanalysis') for which this peak definition applies
%   - subjname: subject for whom this model applies, '*' for all
%   - trialname: name of the trial for which this peak definition applies, '*' for all
%   - direction: specificying wether indicating whether we are looking for a peak above ('positive') or below ('negative') the baseline, 
%     or in any direction ('absolute')
%   - toi: specifies the timewindow in which we are looking for the peak
%       - [begin end] in seconds
%       - pre-specified string: 'all', 'prestim', 'poststim'
%   - peakname: string label for peak 
%       MUST be unique within- and across-trials!
%       MUST NOT contain whitespace, underscore, dash or other character not valid for variable name!
%
% Examples
%aap=aas_add_meeg_peakdef(aap,'aamod_meeg_timelockanalysis','*','STIM','positive',[0.3 0.5],'P300')
%aap=aas_add_meeg_peakdef(aap,'aamod_meeg_timelockanalysis','*','FBNEG','negative','poststim','FRN')

function aap = aas_add_meeg_peakdef(aap, modulename, subjname, trialname, direction, toi, peakname)

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

if any(peakname == '_'), aas_log(aap,true,sprintf('ERROR: Peakname %s has "_"!',peakname)); end

if ~iscell(subjname), subjname = {subjname}; end
if subjname{1} == '*'
    subjname = {aap.acq_details.subjects.subjname};
end

for subj = 1:numel(subjname)
    for m = 1 : length(moduleindex)
        mInd = moduleindex(m);
        
        % find trial
        whichtrialmodel = strcmp({aap.tasksettings.(modulename)(mInd).trialmodel.subject},subjname{subj});
        if ~any(whichtrialmodel), aas_log(aap,true,['ERROR: No trial definition has found for subject: ' subjname{subj}]), end
        whichtrial = strcmp({aap.tasksettings.(modulename)(mInd).trialmodel(whichtrialmodel).model.name},trialname);
        if ~any(whichtrial)
            if strcmp(trialname,'*'),  trialname = {aap.tasksettings.(modulename)(mInd).trialmodel(whichtrialmodel).model.name};
            else, aas_log(aap,true,['ERROR: No trial definition has found: ' trialname]); end
        end
        
        % find model that corresponds and add peakdef to this if it exists
        whichpeakdef = strcmp({aap.tasksettings.(modulename)(mInd).peakdef.subject},subjname{subj}) & ...
            strcmp({aap.tasksettings.(modulename)(mInd).peakdef.trial},trialname);
        
        if (~any(whichpeakdef))
            emptydef = aap.tasksettings.(modulename)(mInd).peakdef(1); % The first one is usually empty, makes for a good template in case the structure changes
            emptydef.subject = subjname{subj};
            emptydef.trial = trialname;
            emptydef.peakdef.direction = direction;
            emptydef.peakdef.toi = toi;
            emptydef.peakdef.name = peakname;
            aap.tasksettings.(modulename)(mInd).peakdef(end+1) = emptydef;
        else
            aap.tasksettings.(modulename)(mInd).peakdef(whichpeakdef).peakdef(end+1).direction=direction;
            aap.tasksettings.(modulename)(mInd).peakdef(whichpeakdef).peakdef(end).toi=toi;
            aap.tasksettings.(modulename)(mInd).peakdef(whichpeakdef).peakdef(end).name=peakname;
        end
    end
end
