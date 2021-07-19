% Adds an contrast for a model
% FORMAT function aap = aas_addcontrast(aap, modulename, subjname, format, vector, conname, contype)
%   - aap: aap structure with parameters and tasklist
%   - modulename: name of module (e.g.,'aamod_firstlevel_contrasts') for which this contrast applies
%   - subjname: subject for whom this model applies
%   - format: format for specificying session(s) for which this contrast applies
%       - "sameforallsessions" -> vector contains contrast to be applied to all sessions (within-subject contrast)
%       - "singlesession:<session name>" -> vector contains contrast for just one session
%       - "sessions:<session name>[+<session name>[...]]" -> vector contains contrast for some sessions
%       - "sessions:<weight>x<session name>[|<weight>x<session name>[...]]" -> vector contains contrast for some (differently) weighted sessions
%           N.B.: You have to use session names with UPPERCASE letters only!
%       - "uniquebysession" -> within-subject contrast that separately specifies contrast for every session
%   - vector: contrast
%       - vector containing the weights for each regressor (of interest)
%       - string defining the weights and the regressors (only for "singlesession:<session name>" and "sameforallsessions")
%           format: <weight>x<regressor name>[<main ('m') or parametric ('p')><number of basis/parametric function>] (e.g. '+1xENC_DISPL|-1xENC_FIXAT' or '+1xENC_DISPLm1|-1xENC_FIXATp3')
%           N.B.: You have to use regressor names with UPPERCASE letters only!
%           number of basis/parametric functions can be, e.g.:
%             - "m2": 2nd order of the basis function (e.g. temporal derivative of canocical hrf) of the main regressor
%             - "p3": depending on the number of basis functions and the order of expansions of parametric modulators (they are arranged after each other)
%                     - dispersion derivative of the 1st order polynomial expansion of the first parametric modulator
%                     - 2nd order polynomial expansion of the first parametric modulator
%                     - 3rd order polynomial expansion of the first parametric modulator
%                     - 1st order polynomial expansion of the second parametric modulator
%                     - 2nd order polynomial expansion of the second parametric modulator
%                     - 1st order polynomial expansion of the third parametric modulator
%   - conname: string label for contrast 
%       Must be unique within- and across-sessions!
%       If empty, then aamod_firstlevel_contrasts will create one.
%   - contype="T" or "F" (defaults to "T")
%   - automatic_movesandmeans=1 or 0, add means & moves to contrast automatically?
%
% Examples
%aap=aas_addcontrast(aap,'aamod_firstlevel_contrasts','*','singlesession:avtask',[0 0 1 1 1])
%aap=aas_addcontrast(aap,'aamod_firstlevel_contrasts','*','singlesession:avtask',[0 0 1 1 1],'stimulus-fixation','T')

function aap = aas_addcontrast(aap, modulename, subjname, format, vector, conname, contype)

% Regexp for number at the end of a module name, if present in format _%05d (e.g. _00001)
m1 = regexp(modulename, '_\d{5,5}$');

% Or, we could use '_*' at the end of the module name to specify all modules with that name
m2 = regexp(modulename, '_\*$');

% Or, we might specify certain modules with  '_X/X/X' (e.g. _00001/00002/00004)
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

if (~exist('conname','var'))
    conname=[];
end;
if (~exist('contype','var') || isempty(contype))
    contype='T';
end;

[format, rem]=strtok(format,':');
switch format
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

% warning
if ischar(vector)
    aas_log(aap,false,'WARNING: You specified the contrast with regressor names.');
    aas_log(aap,false,'    Make sure that you use regressor names with UPPERCASE letters only!',aap.gui_controls.colours.warning);
end

if ~iscell(subjname), subjname = {subjname}; end
if subjname{1} == '*'
    subjname = {aap.acq_details.subjects.subjname};
end

for subj = 1:numel(subjname)
    % check if (any of) the session(s) of the subject exist
    if ~strcmp(format,'uniquebysession')
        sessnames = session.names;
        if isempty(sessnames), sessnames = {aap.acq_details.sessions.name}; end % sameforallsessions
        havesess = true(1,numel(sessnames));
        for s = 1:numel(sessnames)
            sess = find(strcmp({aap.acq_details.sessions.name},sessnames{s}));
            if isempty(sess)
                % avoid cryptic crashes in aas_get_series
                aas_log(aap,true,sprintf(...
                    'did not find sessname %s in {aap.acq_details.sessions.name}',...
                    sessnames{s}));
            end
            [junk, mriser] = aas_get_series(aap,'functional',subj,sess);
            if isempty(mriser) || (isnumeric(mriser) && ~mriser), havesess(s) = false; end
        end
        if ~any(havesess), continue; end
    end
    
    % find model that corresponds and add contrast to this if it exists
    for m = 1 : length(moduleindex)
        
        mInd = moduleindex(m);
        
        whichcontrast=strcmp({aap.tasksettings.(modulename)(mInd).contrasts.subject},subjname{subj});
        if (~any(whichcontrast))
            emptycon=aap.tasksettings.(modulename)(mInd).contrasts(1); % The first one is usually empty, makes for a good template in case the structure changes
            emptycon.subject=subjname{subj};
            emptycon.con.format=format;
            emptycon.con.vector=vector;
            emptycon.con.session=session;
            emptycon.con.type=contype;
            emptycon.con.name=conname;
            aap.tasksettings.(modulename)(mInd).contrasts(end+1)=emptycon;
        else
            
            aap.tasksettings.(modulename)(mInd).contrasts(whichcontrast).con(end+1).format=format;
            aap.tasksettings.(modulename)(mInd).contrasts(whichcontrast).con(end).vector=vector;
            aap.tasksettings.(modulename)(mInd).contrasts(whichcontrast).con(end).session=session;
            aap.tasksettings.(modulename)(mInd).contrasts(whichcontrast).con(end).type=contype;
            aap.tasksettings.(modulename)(mInd).contrasts(whichcontrast).con(end).name=conname;
        end
    end
end
