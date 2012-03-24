function aas_meg_ReportCurrentJob(aap,task,subblock)

% task
tasks=char(fieldnames(aap.MEG.TaskSettings));
dots=repmat('.',[1 size(tasks,2)-length(task)]);
txt=sprintf('\nRunning:%s%s',dots,task);

% subject
if length(subblock)>0
    subs=aap.acq_details.subjects;
    try dots=repmat('.',[1 size(char(subs(1,:,:)),2)-length(aap.acq_details.subjects{subblock{1}})]);
    catch dots='...';
    end
    txt=[txt sprintf('; subject:%s%s',dots,aap.acq_details.subjects{subblock{1}})];
else
    txt=[txt '...'];
end

% block
if length(subblock)==2
    blocks=aap.acq_details.sessions;
    dots=repmat('.',[1 size(char(blocks{1,:,:}),2)-length(aap.acq_details.sessions{subblock{2}})]);
    txt=[txt sprintf('; block:%s%s',dots,aap.acq_details.sessions{subblock{2}})];
else
    txt=[txt '...'];
end

fprintf(txt);
return
