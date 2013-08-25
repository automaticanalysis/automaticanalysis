function aas_recursivekill(aap,hostname,pid)

% Construct a list of all of my processes and their parents
[s w]=aas_shell(['ssh -x ' hostname ' "ps -ef | grep `whoami`"']);
pidlist=[];
while (~isempty(w))
    [nextline w]=strtok(w,[13 10]);
    [field rem]=strtok(nextline);
    [mypid rem]=strtok(rem);
    [myppid rem]=strtok(rem);
    pidlist=[pidlist;str2num(mypid),str2num(myppid)];
end;


% Recursively search out all of the child processes
childlist=[];
for i=1:length(pid)
    childlist=union(childlist,findchildren(pidlist,pid(i),pid(i)));
end;
childlist=flipud(childlist); % kill children first

aas_log(aap,false,sprintf('PARALLEL: terminating processes on %s with pids %s', hostname, sprintf('%d\t',childlist)));

% Now kill all of the processes
for i=1:length(childlist)
    cmd=['ssh -x ' hostname ' "kill -9 ' num2str(childlist(i)) '"'];
    [s w]=unix(cmd);
end;

end


% Recursive search for children
function [childlist]=findchildren(allpids, pid, childlist)
mychildren=allpids(find(allpids(:,2)==pid),1);
childlist=union(childlist,mychildren);
for i=1:length(mychildren)
    childlist=findchildren(allpids,mychildren(i),childlist);
end;
end