function k = aas_getsourcestage(aap,sourcemod,stream)

currmod = '';
currmodnum = aap.tasklist.currenttask.modulenumber;

if nargin == 2 % based on branch
    currbranch = aap.tasklist.currenttask.extraparameters.aap.directory_conventions.analysisid_suffix;
    for currmodnum = find(strcmp({aap.tasklist.main.module.name}, sourcemod))
        srcbranch = aap.tasklist.main.module(currmodnum).extraparameters.aap.directory_conventions.analysisid_suffix;
        if strcmp(currbranch,srcbranch), break; end
    end
    if ~strcmp(currbranch,srcbranch), currmodnum = aap.tasklist.currenttask.modulenumber; end

else % based on stream    
    while ~strcmp(currmod,sourcemod)
        inpstreams = aap.internal.inputstreamsources{currmodnum}.stream;
        try
            currmodnum = inpstreams(strcmp({inpstreams.name},stream)).sourcenumber;
        catch
            aas_log(aap,true,sprintf('stream "%s" cannot be tracked beyond module "%s"',stream,currmod))
        end
        currmod = aap.tasklist.main.module(currmodnum).name;
    end
end

k = currmodnum;
