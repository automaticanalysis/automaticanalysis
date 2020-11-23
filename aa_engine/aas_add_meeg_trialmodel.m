% Adds a model for trials
% FORMAT function aap = aas_add_meeg_trialmodel(aap, modulename, subjname, sessspec, eventspec, samplevector, modelname)
%   - aap: aap structure with parameters and tasklist
%   - modulename: name of module (e.g.,'aamod_meeg_timelockanalysis') for which this trialmodel applies
%   - subjname: subject for whom this model applies, '*' for all
%   - sessspec: specificying session(s) for which this trialmodel applies
%       - "sameforallsessions" -> vector contains trialmodel to be applied to all sessions (within-subject averaging)
%       - "singlesession:<session name>" -> vector considers just one session
%       - "sessions:<session name>[+<session name>[...]]" -> vector contains trialmodel for some sessions
%       - "sessions:<weight>x<session name>[|<weight>x<session name>[...]]" -> vector contains trialmodel for some (differently) weighted sessions
%           N.B.: You have to use session names with UPPERCASE letters only!
%   - eventspec: specificying trial(s) for which this trialmodel applies
%       - vector containing the weights for each trialtype
%       - string defining the wights and the regressors
%           format: <weight>x<regressor name> (e.g. '+1xFB_POS|-1xFB_NEG')
%           N.B.: You have to use event names with UPPERCASE letters only!
%   - samplevector containing the weights for each sample points
%       - function handle to a function converting sample point to weight. E.g. @(x) ones(1,numel(x)) for an average (i.e. weight = 1 for all sample)
%       - pre-specified string
%           - 'avg' - average
%           - 'cont' - contrinuous
%           - 'segmentavg' - average within segments and turn averages into continuous 
%   - modelname: string label for trialmodel 
%       MUST be unique within- and across-sessions!
%       MUST NOT contain whitespace, underscore, dash or other character not valid for variable name!
%
% Examples
%aap=aas_add_meeg_trialmodel(aap,'aamod_meeg_timelockanalysis','*','singlesession:run1',+1xSTIM',@(x) x,'STIMincrease')
%aap=aas_add_meeg_trialmodel(aap,'aamod_meeg_timelockanalysis','*','singlesession:run1','+1xFBPOS|-1xFBNEG','avg','FBPOSminusNEG')

function aap = aas_add_meeg_trialmodel(aap, modulename, subjname, sessspec, eventspec, samplevector, modelname)

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

if any(modelname == '_'), aas_log(aap,true,sprintf('ERROR: Modelname %s has "_"!',modelname)); end

[sessspec, rem]=strtok(sessspec,':');
switch sessspec
    case {'singlesession','sessions'}
        sessstr = strtok(rem,':');
        if any(sessstr == '|') % weighted
            aas_log(aap,false,'WARNING: You specified weights for sessions.');
            aas_log(aap,false,'    Make sure that you use session names with UPPERCASE letters only!',aap.gui_controls.colours.warning);
            cs = textscan(sessstr,'%fx%s','Delimiter','|');
            session.names = cs{2};
            session.weights = cs{1};
        else % simple
            cs = textscan(strtok(rem,':'),'%s','delimiter','+'); 
            session.names = cs{1};
            session.weights = ones(1,numel(session.names));
        end
    otherwise
        session.names=[];
        session.weights = [];
end

if ischar(eventspec)
    if any(eventspec == 'x') % weighted
        aas_log(aap,false,'WARNING: You specified weights for events.');
        aas_log(aap,false,'    Make sure that you use event names with UPPERCASE letters only!',aap.gui_controls.colours.warning);
        cs = textscan(eventspec,'%fx%s','Delimiter','|');
        eventspec = [];
        eventspec.names = cs{2};
        eventspec.weights = cs{1};
    else % simple
        cs = textscan(strtok(eventspec,':'),'%s','delimiter','+');
        eventspec = [];
        eventspec.names = cs{1};
        eventspec.weights = ones(1,numel(eventspec.names));
    end
end

if ~iscell(subjname), subjname = {subjname}; end
if subjname{1} == '*'
    subjname = {aap.acq_details.subjects.subjname};
end

for subj = 1:numel(subjname)
    % check if (any of) the session(s) of the subject exist
    sessnames = session.names;
    if isempty(sessnames), sessnames = {aap.acq_details.meeg_sessions.name}; end % sameforallsessions
    havesess = true(1,numel(sessnames));
    for s = 1:numel(sessnames)
        sess = find(strcmp({aap.acq_details.meeg_sessions.name},sessnames{s}));
        if isempty(sess)
            % avoid cryptic crashes in aas_get_series
            aas_log(aap,true,sprintf(...
                'did not find sessname %s in {aap.acq_details.meeg_sessions.name}',...
                sessnames{s}));
        end
        [junk, meegser] = aas_get_series(aap,'MEEG',subj,sess);
        if isempty(meegser) || (isnumeric(meegser) && ~meegser), havesess(s) = false; end
    end
    if ~any(havesess), continue; end
    
    % find model that corresponds and add trialmodel to this if it exists
    for m = 1 : length(moduleindex)
        
        mInd = moduleindex(m);
        
        whichmodule = strcmp({aap.tasksettings.(modulename)(mInd).trialmodel.subject},subjname{subj});
        if (~any(whichmodule))
            emptymod = aap.tasksettings.(modulename)(mInd).trialmodel(1); % The first one is usually empty, makes for a good template in case the structure changes
            emptymod.subject=subjname{subj};
            emptymod.model.samplevector=samplevector;
            emptymod.model.session=session;
            emptymod.model.event=eventspec;
            emptymod.model.name=modelname;
            aap.tasksettings.(modulename)(mInd).trialmodel(end+1)=emptymod;
        else
            
            aap.tasksettings.(modulename)(mInd).trialmodel(whichmodule).model(end+1).samplevector=samplevector;
            aap.tasksettings.(modulename)(mInd).trialmodel(whichmodule).model(end).session=session;
            aap.tasksettings.(modulename)(mInd).trialmodel(whichmodule).model(end).event=eventspec;
            aap.tasksettings.(modulename)(mInd).trialmodel(whichmodule).model(end).name=modelname;
        end
    end
end
