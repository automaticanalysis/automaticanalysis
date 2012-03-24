% AA initialisation module - check aap has been correctly set up by user
% Uses a copy of the parameters made taken straight after the recipe
%  construction in aa_init
% If aap structure is different, mistyping has probably occurred and an
% error will be generated
% Rhodri Cusack MRC CBU Cambridge 2004


function [aap,resp]=aamod_checkparameters(aap,task)

resp='';

switch task
    case 'domain'
        resp='study';   % this module needs to be run once per study
    case 'description'
        resp='Check user parameters';
    case 'summary'
        resp=['Check no parameter names mistyped'\n'];
    case 'doit'
        
        aap_afteruser=aap;
        
        % mask the aap_beforeuserchange structure
        aap_beforeuser=aap.aap_beforeuserchanges;
        aap_afteruser.aap_beforeuserchanges=[];
        
        % mask the internal fields which are allowed to change
        aap_beforeuser.internal=[];
        aap_afteruser.internal=[];
        
        aas_recursivecompare(aap,aap_beforeuser,aap_afteruser,'aap');
        
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;    

return;


% RECURSIVE COMPARISON OF STRUCTURES BEFORE AND AFTER USER CHANGES

function aas_recursivecompare(aap,aap_beforeuser,aap_afteruser,pth);

bef_fieldnames=fieldnames(aap_beforeuser);
for i=1:length(bef_fieldnames)
    newpth=[pth '.' bef_fieldnames{i}];
    if (~isfield(aap_afteruser,bef_fieldnames{i}))
        aas_log(aap,1,sprintf('Field "%s" is present in recipe before user changes but was not after\n',newpth));
    end;
    if (isstruct(getfield(aap_afteruser,bef_fieldnames{i})) & ~isstruct(getfield(aap_beforeuser,bef_fieldnames{i})))
        aas_log(aap,1,sprintf('Field "%s" is a structure after user changes but was not before',newpth))
    end;
    if (~isstruct(getfield(aap_afteruser,bef_fieldnames{i})) & isstruct(getfield(aap_beforeuser,bef_fieldnames{i})))
        aas_log(aap,1,sprintf('Field "%s" was a structure before user changes but was not after',newpth))
    end;
end;

aft_fieldnames=fieldnames(aap_afteruser);
for i=1:length(aft_fieldnames)
    newpth=[pth '.' aft_fieldnames{i}];
    if (~isfield(aap_beforeuser,aft_fieldnames{i}))
        aas_log(aap,1,sprintf('Field "%s" is present in recipe after user changes but not before\n',newpth));
    end;
    if (isstruct(getfield(aap_afteruser,aft_fieldnames{i})) & ~isstruct(getfield(aap_beforeuser,aft_fieldnames{i})))
        aas_log(aap,1,sprintf('Field "%s" is a structure after user changes but was not before',newpth))
    end;
    if (~isstruct(getfield(aap_afteruser,aft_fieldnames{i})) & isstruct(getfield(aap_beforeuser,aft_fieldnames{i})))
        aas_log(aap,1,sprintf('Field "%s" was a structure after user changes but was not before',newpth))
    end;
    
    if (isstruct(getfield(aap_afteruser,aft_fieldnames{i})))
        aas_recursivecompare(aap,getfield(aap_beforeuser,aft_fieldnames{i}),getfield(aap_afteruser,aft_fieldnames{i}),newpth);
    end;
end;


