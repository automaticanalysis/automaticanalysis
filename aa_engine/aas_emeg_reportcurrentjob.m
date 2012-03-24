function aas_emeg_reportcurrentjob(aap,task,subblock,action)

if nargin<4; action='Running'; end;

if isnumeric(subblock); subblock={subblock}; end % allow variable input types

% task
tasks=char(fieldnames(aap.tasksettings));
dots=repmat('.',[1 size(tasks,2)-length(task)]);
txt=sprintf('\n%s:%s%s',action,dots,task);

% subject
if length(subblock)>0
    try
        if strcmp('internal',aap.schema.tasksettings.(task).ATTRIBUTE.domain)
            txt=[txt sprintf(' #%g\n',subblock{1})]; % for any internal domain
        else
            subs=struct2cell(aap.acq_details.subjects);
            subname=aap.acq_details.subjects(subblock{1}).megname;
            if isempty(subname);
                subname=aap.acq_details.subjects(subblock{1}).mriname;
            end
            try dots=repmat('.',[1 size(char(subs(1,:,:)),2)-length(subname)]);
            catch dots='...';
            end
            txt=[txt sprintf('; subject:%s%s\n',dots,aap.acq_details.subjects(subblock{1}).megname)];
        end
    catch
        % if manually running a job on an aap structure that didn't include it then can ignore this
    end
else
    txt=[txt '...'];
end

% block
if length(subblock)==2
    blocks=struct2cell(aap.acq_details.sessions);
    dots=repmat('.',[1 size(char(blocks(1,:,:)),2)-length(aap.acq_details.sessions(subblock{2}).name)]);
    txt=[txt sprintf('; block:%s%s\n',dots,aap.acq_details.sessions(subblock{2}).name)];
else
    txt=[txt '...'];
end

fprintf(txt);
return