function aap=aas_add_meg_event(aap,modulename,subject,session,eventname,eventvalue,trialshift)

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

subject = aas_mriname2subjname(subject);

event.conditionlabel=eventname;
event.eventvalue=eventvalue;
event.trialshift=trialshift;

% find models that corresponds and add events if they exist
for m = 1 : length(moduleindex)
    
    mInd = moduleindex(m);
    
    whichmodel=[strcmp({aap.tasksettings.(modulename)(mInd).condition.subject},subject)] & [strcmp({aap.tasksettings.(modulename)(mInd).condition.session},session)];
    if (~any(whichmodel))
        emptymod=[];
        emptymod.subject=subject;
        emptymod.session=session;
        emptymod.event = event;
        emptymod.event.eventtype = aap.tasksettings.(modulename)(mInd).eventtype;
        aap.tasksettings.(modulename)(mInd).condition(end+1)=emptymod;
    else
        event.eventtype = aap.tasksettings.(modulename)(mInd).eventtype;
        aap.tasksettings.(modulename)(mInd).condition(whichmodel).event(end+1) = event;
    end
    
end