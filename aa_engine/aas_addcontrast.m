% Adds an contrast for a model
% function aap=aas_addcontrast(aap,modulename,subject,format,vector,conname,contype,automatic_movesandmeans)
%
% modulename=name of module (e.g.,'aamod_firstlevel_contrasts') for which this
%   contrast applies
% subject=subject for whom this model applies
% format=format for contrast specification, one of:
%  "sameforallsessions" - vector contains contrast to be applied to all sessions
%  "singlesession:[sessionname]" - vector contains contrast for just one
%     session, all other sessions will be set to 0. [sessionname] should be
%     replaced with name of that session.
%  "uniquebysession" - long contrast string that separately specifies
%     contrast for every session
% vector=vecor containing the contrast
% conname=cell string label for contrast (if empty, then aamod_firstlevel_contrast will create one)
% contype="T" or "F" (defaults to "T")
% automatic_movesandmeans=1 or 0, add means & moves to contrast
%   automatically?
%
% Examples
%aap=aas_addcontrast(aap,'aamod_firstlevel_contrasts','*','singlesession:avtask',[0 0 1 1 1])
%aap=aas_addcontrast(aap,'aamod_firstlevel_contrasts','*','singlesession:avtask',[0 0 1 1 1],'stimulus-fixation','T')

function aap=aas_addcontrast(aap,modulename,subject,format,vector,conname,contype,automatic_movesandmeans)

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

% % Get number from end of module name if present in format _%05d (e.g, _00001)
% if (length(modulename>6))
%     moduleindex=str2num(modulename(length(modulename)-4:end));
%     if (~strcmp(['_' sprintf('%05d',moduleindex)],modulename(length(modulename)-5:end)))
%         moduleindex=1;
%     else
%         modulename=modulename(1:length(modulename)-6);
%     end;
% else
%     moduleindex=1;
% end;


if (~exist('conname','var'))
    conname=[];
end;

if (~exist('contype','var') || isempty(contype))
    contype='T';
end;

if (~exist('automatic_movesandmeans','var') || isempty(automatic_movesandmeans))
    automatic_movesandmeans=true;
end;

[format rem]=strtok(format,':');
if (strcmp(format,'singlesession'))
    session=strtok(rem,':');
else
    session=[];
end;

% find model that corresponds and add contrast to this if it exists
for m = 1 : length(moduleindex)
    
    mInd = moduleindex(m);
    
    whichcontrast=[strcmp({aap.tasksettings.(modulename)(mInd).contrasts.subject},subject)];
    if (~any(whichcontrast))
        emptycon=[];
        emptycon.subject=subject;
        emptycon.automatic_movesandmeans=automatic_movesandmeans;
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

