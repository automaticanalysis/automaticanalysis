function [aap]=aas_clear_jobq(aap,username,wheretoprocess)
if (~exist('wheretoprocess','var') || isempty(wheretoprocess))
    wheretoprocess=aap.options.wheretoprocess;
end;
switch (wheretoprocess)
    case 'localsingle'
        aas_log(aap,false,'aas_clear_jobq - no queue as mode localsingle');
    case 'aws'
        aap=aws_setupqnames(aap,username);
        q=aaq_aws(aap);
        q.clear_jobq();
    case 'localparallel'
        aas_log(aap,false,'aas_clear_jobq - not yet implemented for localparallel');
end;
end