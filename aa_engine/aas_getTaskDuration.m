function elapsed = aas_getTaskDuration(Task)
try
    DisplayHelper = parallel.internal.display.DisplayHelper('Task');
catch E
    if strcmp(E.identifier,'MATLAB:invalidType'), DisplayHelper = parallel.internal.display.DisplayHelper(4);
    else rethrow(E); end
end

if isprop(Task,'StartDateTime')
    startTime = Task.StartDateTime;
    finishTime = Task.FinishDateTime;
elseif isprop(Task,'StartTime')
    startTime = Task.StartTime;
    finishTime = Task.FinishTime;
else
    aas_log(obj.aap,true,'Time-related property of Task class not found!')
end

try
    elapsed = DisplayHelper.getRunningDuration(startTime,finishTime);
catch E
    if strcmp(E.identifier,'MATLAB:invalidType')
        elapsed = DisplayHelper.getRunningDuration(datetime2java(startTime),datetime2java(finishTime));
    else rethrow(E);
    end
end

end

function str = datetime2java(dt)
    str = [datestr(dt,'ddd mmm dd HH:MM:SS ') dt.TimeZone datestr(dt,' yyyy')];
end