function aap = aas_addcontrast(aap, modulename, subjname, sessionspec, condef, conname, contype)
%
% Add a contrast to a model
%
% inputs
%
%   aap - aap structure
%
%   modulename - module(s) for which this contrast applies. Can be a
%   single module, or multiple modules identified using A regex string.
%   See examples below.
%
%   subjname - subjects for which this contrast applies. Can be a single
%   subject ID (string), a cell array of multiple subject IDs, or a wildcard
%   ("*") to apply constrast to all subjets. See examples below.
%
%   sessionspec - sessions for which this contrast applies. 
%   The following options are recognized:
%
%   "sameforallsessions" 
%       - contrast applies to all sessions defined in the model
%
%   "singlesession:<session name>" 
%       - contrast applies to the single named session
%
%   "sessions:<session name>[+<session name>[...]]" 
%       - contrast applies to multiple named sessions
%
%   "sessions:<weight>x<session name>[|<weight>x<session name>[...]]" 
%       - contrast apples to multiple sessions differently weighted
%         N.B.: this option requires session names are UPPERCASE
%
%   "uniquebysession" 
%       - contrast is appled across *all* sessions
%
%   See examples below
%
%   "condef" - contrast definition. Two options are recognized:
%
%   1) (numeric) vector. This will be zero padded to the number of columns
%   in each session, or to the number of columns in the design matrix
%   for "uniquebysession"
%
%   2) string specifier in the form 
%
%       <weight>x<REGRESSORNAME>[|<weight>x<REGRESSORNAME>...]
%
%   where "weight" is a signed string expression (e.g., +1, -1, etc)
%
%   this can optionally include main or parametric order:
%
%       <weight>x<REGRESSORNAMEmN>
%       <weight>x<REGRESSORNAMEpN>
%
%   that is:
%
%       <weight>x<regressor name>[<main ('m') or parametric ('p')><number of basis/parametric function>]
%
%   for example:
%   
%       '+1xENC_DISPLm1|-1xENC_FIXATp3')
%
%   where the number of basis/parametric functions can be, e.g.:
%
%     - "m2": 2nd order of the basis function (e.g. temporal derivative of canocical hrf) of the main regressor
%     - "p3": depending on the number of basis functions and the order of expansions of parametric modulators 
%         - dispersion derivative of the 1st order polynomial expansion of the first parametric modulator
%         - 2nd order polynomial expansion of the first parametric modulator
%         - 3rd order polynomial expansion of the first parametric modulator
%         - 1st order polynomial expansion of the second parametric modulator
%         - 2nd order polynomial expansion of the second parametric modulator
%         - 1st order polynomial expansion of the third parametric modulator
%
%   N.B.: use of this option requires event names to be specified in UPPERCASE
%
%   Additional examples are provided below.
%
%   conname - contrast name. Must be unique within and across sessions.
%   If empty, then aamod_firstlevel_contrasts will create a name.
%
%   contype = "T" or "F" (defaults to "T")
%
% Notes
%
% 1) call aas_addcontrast once for each contrast to be defined in the model
%
% 2) avoid the use of spaces or special characters (~!@#$%^&*()_+}{|: etc)
%    in contrast and regressor names. Specifically, avoid the use of > and
%    > in the contrast name. Addtionally, while only a few options require
%    session and/or regressor names to be in uppercase, it is generally
%    good programming practice to always use uppercase for session and
%    regressor names
%
% 3) an easy mistake is to forget to include the aap struct as a return
%    variable when calling addcontrast in your userscript. That is,
%    be sure to do: aap = aas_addcontrast(aap,.... and not
%    aas_addcontrast(aap,...).
%
% 4) note the standard first level contrast module in aa is named
%    aamod_firstlevel_contrasts (plural), not aamod_firstlevel_contrast
%
% 5) As a rule, differential contrasts should sum to zero to avoid
%    scaline issues. See any good text on the GLM for details.
%
% ------------------------------------------------------------------------------------------------------------------------------------
%
% Examples
%
% Example-1: Assume the model contains two sessions named SESS01 and SESS02, each with two regressors of interest 
% called R1 and R2 and the six motion parameters (x,y,x,r,p,j) as nuisance regressors.
%
% 1a) This call defines an omnibus contrast for R1>0 for all subjects:
% 
%    aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts', '*', 'sameforallsessions', [1], 'R1_G_0', 'T');
%
%   Matlab might warn [1] is unncecessary and can be replaced by 1:
%
%    aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts', '*', 'sameforallsessions', 1, 'R1_G_0', 'T');
%
% Note in both examples, we avoided the use of ">" in the contrast name, instead using "_G_"
%
% 1b) This call defines the differential contrast R1>R2:
%
%    aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts', '*', 'sameforallsessions', [1 -1], 'R1_G_R2', 'T');
%
%   We can define this contrast for sub-01 only:
%
%    aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts', {'sub-01'}, 'sameforallsessions', [1 -1], 'R1_G_R2', 'T');
%
%   We can define this contrast for sub-01 and SESS01 only
%
%    aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts', {'sub-01'}, 'singlesession:SESS01', [1 -1], 'R1_G_R2', 'T');
%
%   If the tasklist contains more than one instance of aamod_firstlevel_contrasts, we can define this contrast for all
%   instances using a wildcard:
%
%    aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts_*', {'sub-01'}, 'session:SESS01', [1 -1], 'R1_G_R2', 'T');
%
%   or, we could define this contrast only for a specific instance by including the module numerical suffix:
%
%    aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts_00002', {'sub-01'}, 'session:SESS01', [1 -1], 'R1_G_R2', 'T');
%
%   would define the contrast only for the second instance of aamod_firstlevel_contrast. Note the module suffix
%   is FIVE characters (i.e. _00001, not _0001 or _001 or _01 or _1)
%
% 1c) We could define the differential constast R1>R2 (1a) using uniquebysession:
%
%   cvec = [1 -1   0 0 0 0 0 0   1 -1   0 0 0 0 0 0    0 0 ];
%   aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts', '*', 'uniquebysession', cvec, 'R1_G_R2', 'T');
%
%   note the contrast vector must include entries for all columns of the design matrix when using uniquebysession.
%   This includes entries for the six motion parameters in each session, and two final zeros for the two constant 
%   terms automatically added to the model.
%
%   However, aa will automatically pad the lenght of the contrast vector if necessary. So we could define: 
%
%   cvec = [1 -1   0 0 0 0 0 0   1 -1 ];
%   aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts', '*', 'uniquebysession', cvec, 'R1_G_R2', 'T');
%
%   However, it's usually easier to define a uniquebysession contrast using a string specifier:
%
%   aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts', '*', 'uniquebysession', '+1xR1|-1xR2', 'R1_G_R2', 'T');
%
% This is of course a silly example, because it's much simpler to define R1_G_R2 using sameforallsessions as in 1b. 
% The usefulness of uniquebysession becomes apparent when different sessions contain DIFFERENT regressors. See
% next example.
%
% Example-2 Consider a model with two sessions: SESS01 has events LC (listen to a Clear word) and LN (Listen to a word 
% presented in Noise). SESS02 has events RC (Repeat a word presented in Clear) and RN (Repeat a word presented in Noise)
%
% We could define the contrast LC>LN and RC>RN using the single session specifier:
%
% aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts', '*', 'session:SESS01', [1 -1], 'LC_G_LN', 'T');
% aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts', '*', 'session:SESS02', [1 -1], 'RC_G_RN', 'T');
%
% However, we can't use this approach to define, say, LC>RC, because the contrast vector spans 
% sessions and the columns mean different things in different sessions. The following is WRONG:
%
% aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts', '*', 'sameforallsessions', [1 -1], 'LC_G_RC', 'T'); *** WRONG
%
% Instead we must use uniquebysession:
%
%   cvec = [1 0   0 0 0 0 0 0   -1 0    0 0 0 0 0 0    0 0 ];
%   aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts', '*', 'uniquebysession', cvec, 'LC_G_RC', 'T');
%   cvec = [0 1   0 0 0 0 0 0   0 -1    0 0 0 0 0 0    0 0 ];
%   aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts', '*', 'uniquebysession', cvec, 'LN_G_RN', 'T');
%
% (As before, we assume the six motion parameters are included in the model as nuisance regressors.)
%
% However, a better approach is to define these contrasts using a string specifier:
%
%   aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts', '*', 'uniquebysession', '+1xLC|-1xRC', 'LC_G_RC', 'T');
%   aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts', '*', 'uniquebysession', '+1xLN|-1xRN', 'LN_G_RN', 'T');
%
% Footnote: This is essential when including frame censoring in your pipeline, because (unlike the six motion paramters) 
% you don't know a priori how many nuisance regressors will appear in each session. As such, you would not be able to
% define a contrast spanning multiple sessions using a zero-padded vector.
%
% Remember, you cannot use the string specifier option unless your regressor names are UPPERCASE
%
% ------------------------------------------------------------------------------------------------------------------------------------

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
end
if (~exist('contype','var') || isempty(contype))
    contype='T';
end

[sessionspec, rem]=strtok(sessionspec,':');
switch sessionspec
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

if ischar(condef)
    aas_log(aap,false,'INFO: You specified the contrast with regressor names.');
    aas_log(aap,false,'Regressor names must be UPPERCASE only.',aap.gui_controls.colours.warning);
end

if ~iscell(subjname), subjname = {subjname}; end
if subjname{1} == '*'
    subjname = {aap.acq_details.subjects.subjname};
end

for subj = 1:numel(subjname)
    % check if (any of) the session(s) of the subject exist
    if ~strcmp(sessionspec,'uniquebysession')
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
            [~, mriser] = aas_get_series(aap,'functional',subj,sess);
            if isempty(mriser) || (isnumeric(mriser) && ~mriser), havesess(s) = false; end
        end
        if ~any(havesess), continue; end
    end
    
    % find model that corresponds and add contrast to this if it exists
    for m = 1 : length(moduleindex)
        
        mInd = moduleindex(m);
        
        whichcontrast=strcmp({aap.tasksettings.(modulename)(mInd).contrasts.subject},subjname{subj});
        if (~any(whichcontrast))
            % The first one is usually empty, makes for a good template in case the structure changes
            emptycon=aap.tasksettings.(modulename)(mInd).contrasts(1); 
            emptycon.subject=subjname{subj};
            emptycon.con.format=sessionspec;
            emptycon.con.vector=condef;
            emptycon.con.session=session;
            emptycon.con.type=contype;
            emptycon.con.name=conname;
            aap.tasksettings.(modulename)(mInd).contrasts(end+1)=emptycon;
        else
            
            aap.tasksettings.(modulename)(mInd).contrasts(whichcontrast).con(end+1).format=sessionspec;
            aap.tasksettings.(modulename)(mInd).contrasts(whichcontrast).con(end).vector=condef;
            aap.tasksettings.(modulename)(mInd).contrasts(whichcontrast).con(end).session=session;
            aap.tasksettings.(modulename)(mInd).contrasts(whichcontrast).con(end).type=contype;
            aap.tasksettings.(modulename)(mInd).contrasts(whichcontrast).con(end).name=conname;
        end
    end
end
