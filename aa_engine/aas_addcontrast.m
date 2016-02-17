% Adds an contrast for a model
% FORMAT function aap = aas_addcontrast(aap, modulename, subjname, format, vector, conname, contype)
%   - aap: aap structure with parameters and tasklist
%   - modulename: name of module (e.g.,'aamod_firstlevel_contrasts') for which this contrast applies
%   - subjname: subject for whom this model applies
%   - format: format for specificying session(s) for which this contrast applies
%       - "sameforallsessions" -> vector contains contrast to be applied to all sessions (within-subject contrast)
%       - "singlesession:<session name>" -> vector contains contrast for just one session
%       - "uniquebysession" -> within-subject contrast that separately specifies contrast for every session
%   - vector: contrast
%       - vector containing the weights for each regressor (of interest)
%       - string defining the wights and the regressors in a format <weight>x<regressor name><main ('m') or parametric ('p')><number of basis/parametric function> (e.g. '+1xENC_DISPLm1|-1xENC_FIXATp3')
%           N.B.: number of basis/parametric functions can be, e.g.:
%             - "m2": 2nd order of the basis function (e.g. temporal derivative of canocical hrf) of the main regressor
%             - "p3": depending on the number of basis functions and the order of expansions of parametric modulators (they are arranged after each other)
%                     - dispersion derivative of the 1st order polynomial expansion of the first parametric modulator
%                     - 2nd order polynomial expansion of the first parametric modulator
%                     - 3rd order polynomial expansion of the first parametric modulator
%                     - 1st order polynomial expansion of the second parametric modulator
%                     - 2nd order polynomial expansion of the second parametric modulator
%                     - 1st order polynomial expansion of the third parametric modulator
%   - conname: string label for contrast (if empty, then aamod_firstlevel_contrasts will create one)
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
if (strcmp(format,'singlesession'))
    session=strtok(rem,':');
else
    session=[];
end;

% find model that corresponds and add contrast to this if it exists
for m = 1 : length(moduleindex)
    
    mInd = moduleindex(m);
    
    whichcontrast=strcmp({aap.tasksettings.(modulename)(mInd).contrasts.subject},subjname);
    if (~any(whichcontrast))
        emptycon=aap.tasksettings.(modulename)(mInd).contrasts(1); % The first one is usually empty, makes for a good template in case the structure changes
        emptycon.subject=subjname;
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

