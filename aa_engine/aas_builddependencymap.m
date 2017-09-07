%% BUILD DEPENDENCY MAP
function [aap]=aas_builddependencymap(aap)

aap.internal.dependenton=cell(length(aap.tasklist.main.module),1);
aap.internal.dependencyof=cell(length(aap.tasklist.main.module),1);

for k1=1:length(aap.tasklist.main.module)
    % Default is no dependencies
    completefirst=[];
    
    % first stage has no dependencies
    if (k1>1)
        
        tobecompletedfirst=aap.tasklist.main.module(k1).tobecompletedfirst;
        
        % Defaults to previous stage
        if (isempty(tobecompletedfirst))
            tobecompletedfirst={'[previous]'};
        end;
        if (~iscell(tobecompletedfirst) && ~isempty(tobecompletedfirst))
            tobecompletedfirst={tobecompletedfirst};
        end;
        
        for k0i=1:length(tobecompletedfirst)
            % previous stage, or something elsE?
            if (strcmp(tobecompletedfirst{k0i},'[previous]'))  % empty now same as [previous]
                completefirst(k0i).stage=k1-1;
            elseif (strcmp(tobecompletedfirst{k0i},'[none]'))   %... so none must be explcitly specified
            else
                for k0=1:k1
                    % allow full path of module to be provided
                    stagetag=aas_getstagetag(aap,k0);
                    if strcmp(stagetag,tobecompletedfirst{k0i})
                        break;
                    end;
                    if (k0==k1)
                        aas_log(aap,1,sprintf('%s is only to be executed after %s, but this is not in the task list',aap.tasklist.main.module(k1).name,tobecompletedfirst{k0i}));
                    end;
                end;
                completefirst(k0i).stage=k0;
            end;
            % now find out what domain the done flags need to cover
            % allow full path of module to be provided
            [stagepath stagename]=fileparts(aap.tasklist.main.module(completefirst(k0i).stage).name);
            completefirst(k0i).sourcedomain=aap.schema.tasksettings.(stagename)(aap.tasklist.main.module(completefirst(k0i).stage).index).ATTRIBUTE.domain;
        end;
        
        aap.internal.dependenton{k1}=completefirst;
        
        % now produce indexing the other way, to map
        % aap.internal.dependencyof
        [stagepath stagename]=fileparts(aap.tasklist.main.module(k1).name);
        thisstage=[];
        thisstage.stage=k1;
        thisstage.domain=aap.schema.tasksettings.(stagename)(aap.tasklist.main.module(k1).index).ATTRIBUTE.domain;
        for k0i=1:length(completefirst)
            aap.internal.dependencyof{completefirst(k0i).stage}=[aap.internal.dependencyof{completefirst(k0i).stage} thisstage];
        end;
    end;
    
end;

