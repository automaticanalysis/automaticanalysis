function [aap]=aas_validatepaths(aap)

validatepaths(aap,aap.schema,aap,'aap');


end

function validatepaths(aap,mystruct,values,nme)
if (isfield(mystruct,'ATTRIBUTE'))
    if (~isfield(mystruct.('ATTRIBUTE'),'ui'))
        aas_log(aap,true,sprintf('Definition of parameter %s in XML does not contain ui type\n',nme));
    end;
    ui=mystruct.('ATTRIBUTE').ui;
%     fprintf('Found %s leaf type %s\n',nme,ui);
    switch(ui)
        case {'dir','filename','dir_part'}
            aap=aas_checkpath(aap,values,nme);
    end;
else
    if (~strcmp(nme,'aap.tasksettings'))
        
        fs=fieldnames(mystruct);
        for f=1:length(fs)
            sub=mystruct.(fs{f});
            if (~iscell(sub))
                sub={sub};
            end;
            for ind=1:length(sub)
                if (isstruct(sub{ind}))
                    for structind=1:length(values.(fs{f}))
                        validatepaths(aap,mystruct.(fs{f}),values.(fs{f})(structind),[nme '.' fs{f}]);
                    end;
                else
                    aas_log(aap,true,sprintf('Definition of field %s.%s in XML does not contain any descriptions (e.g., ui)\n',nme,fs{f}));
                end;
            end;
        end;
    end;
end;

end