function [resp]=aas_getworkerdescription(aap,title,summaryworkerlist)
global aaparallel
resp=['<table><tr><td>' num2str(length(summaryworkerlist)) ' ' title '</td></tr>'];
for i=1:length(summaryworkerlist)
    workerid=summaryworkerlist(i);
    resp=[resp '<tr><td>' sprintf('W%d',workerid)];
    Jalloc=[];
    if (~isempty(aaparallel.workerstatus.(sprintf('worker%d',workerid)).allocatedjobs))
        for j=1:length(aaparallel.workerstatus.(sprintf('worker%d',workerid)).allocatedjobs)
            Jalloc=[Jalloc aaparallel.workerstatus.(sprintf('worker%d',workerid)).allocatedjobs{j}.taskqueueposition];
        end;
    end;
    if (~isempty(aaparallel.workerstatus.(sprintf('worker%d',workerid)).runningjob))    
        Jrun=[aaparallel.workerstatus.(sprintf('worker%d',workerid)).runningjob.taskqueueposition];
    else
        Jrun=[];
    end;
    if (~isempty(Jalloc) || ~isempty(Jrun))
        if (~isempty(Jalloc))
            resp=[resp sprintf('-aJ%d',Jalloc)];
        end;
        if (~isempty(Jrun))
            resp=[resp sprintf('-J%d',Jrun)];
        end;
    end;
    resp=[resp '</td></tr>'];
end;
resp=[resp '</table>'];   
