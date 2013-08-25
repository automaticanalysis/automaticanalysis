% Adds an event to a model
% function aap=aas_addevent(aap,modulename,subject,session,eventname,ons,dur,parametric)
%
% modulename=name of module (e.g.,'aamod_firstlevel_model') for which this
%   event applies
% subject=subject for whom this model applies
% session=session for which this applies
% eventname=name of the stimulus or response event
% ons=event onset times (in scans). Does not need to be sorted
% dur=event durations (in scans), either a single element (if all
%   occurrences have the same duration) or in order that corresponds to ons
% parametric = (multiple) parametric modulator with 3 fields
%   parametric(n).name = name of the modulator
%   parametric(n).P = modulator vector itself
%   parametric(n).h = polynomial expansion
%
% Examples
%  aap=aas_addevent(aap,'aamod_firstlevel_model','*','*','AudVid300',ons,dur);
%  aap=aas_addevent(aap,'aamod_firstlevel_model','*','*','AudVid300',ons,dur,parametric);

function aap=aas_addevent(aap,modulename,subject,session,eventname,ons,dur,parametric)

% Regexp for number at the end of a module name, if present in format _%05d (e.g, _00001)
m1 = regexp(modulename, '_\d{5,5}$');

% Or, we could use '_*' at the end of the module name to specify all modules with that name
m2 = regexp(modulename, '_\*$');

% Or, we might specify specific modules with  '_X/X/X'
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

if (~exist('parametric','var'))
    parametric=[];
end

% sort the onsets, and apply same reordering to dur & parametric
% [AVG] - replacd junk by ons, since we *DO* want to sort the onsets
[ons ind]=sort(ons);
if (length(dur)>1)
    dur=dur(ind);
end;
if (~isempty(parametric))
    % [AVG] reorder parametric modulator even if there's more than one!
    for p = 1:length(parametric)
        parametric(p).P = parametric(p).P(ind);
    end
end

% find models that corresponds and add events if they exist
for m = 1 : length(moduleindex)
    
    mInd = moduleindex(m);
    
    whichmodel=[strcmp({aap.tasksettings.(modulename)(mInd).model.subject},subject)] & [strcmp({aap.tasksettings.(modulename)(mInd).model.session},session)];
    if (~any(whichmodel))
        emptymod=[];
        emptymod.subject=subject;
        emptymod.session=session;
        emptymod.event.name=eventname;
        emptymod.event.ons=ons;
        emptymod.event.dur=dur;
        emptymod.event.parametric=parametric;
        aap.tasksettings.(modulename)(mInd).model(end+1)=emptymod;
    else
        aap.tasksettings.(modulename)(mInd).model(whichmodel).event(end+1).name=eventname;
        aap.tasksettings.(modulename)(mInd).model(whichmodel).event(end).ons=ons;
        aap.tasksettings.(modulename)(mInd).model(whichmodel).event(end).dur=dur;
        aap.tasksettings.(modulename)(mInd).model(whichmodel).event(end).parametric=parametric;
    end
    
end