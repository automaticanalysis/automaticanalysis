%%==============================
% Delete done flags
% function aas_delete_doneflag(aap,stage,inpi,inpj)

function aap=aas_delete_doneflag(aap,stage,inpi,inpj)

% allow full path of module to be provided
[stagepath stagename]=fileparts(aap.tasklist.main.module(stage).name);
index=aap.tasklist.main.module(stage).index;
% retrieve description from module
description=aap.schema.tasksettings.(stagename)(index).ATTRIBUTE.desc;
% find out whether this module needs to be executed once per study, subject or session
domain=aap.schema.tasksettings.(stagename)(index).ATTRIBUTE.domain;

if (~exist('inpi','var'))
    boti=1;
    topi=length(aap.acq_details.subjects);
else
    boti=inpi;
    topi=inpi;
end;

if (~exist('inpj','var'))
    botj=1;
    topj=length(aap.acq_details.sessions);
else
    botj=inpj;
    topj=inpj;
end;

switch (domain)
    case 'study'
        doneflag=aas_doneflag_getpath(aap,stage);
        aas_delete_doneflag_bypath(aap,doneflag);
        
    case 'subject'
        
        for i=boti:topi
            doneflag=aas_doneflag_getpath(aap,i,stage);
            aas_delete_doneflag_bypath(aap,doneflag);
            
        end;
        
    case 'session'
        for i=boti:topi
            for j=botj:topj
                doneflag=aas_doneflag_getpath(aap,i,j,stage);
                aas_delete_doneflag_bypath(aap,doneflag);
                
            end;
        end;
end;
