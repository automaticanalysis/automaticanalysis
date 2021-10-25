% Validate that paths specified in aap struct meet format requirements by
% calling aas_checkpath. Called during aa_doprocessing.
%
% [aap]=aas_validatepaths(aap)
function [aap]=aas_validatepaths(aap)

validatepaths(aap,aap.schema,aap,'aap');


end

function validatepaths(aap,mystruct,values,nme)
if (isfield(mystruct,'ATTRIBUTE'))
    if (~isfield(mystruct.('ATTRIBUTE'),'ui'))
        aas_log(aap,true,sprintf('Definition of parameter %s in XML does not contain ui type\n',nme));
    end
    ui=mystruct.('ATTRIBUTE').ui;
    switch(ui)
        case {'dir','filename','dir_part','dir_list'}
            % allow colons in dir_list type
            allowcolons=strcmp('dir_list',ui);

            if isstruct(values) % multiple values
                f = fieldnames(values);
                values = values.(f{1});
                for i = 1:numel(values)
                    aap=aas_checkpath(aap,values{i},nme,[],allowcolons);
                end
            else
                aap=aas_checkpath(aap,values,nme,[],allowcolons);
            end
        otherwise
            % no checking
    end
else
    if (~strcmp(nme,'aap.tasksettings'))

        fs=fieldnames(mystruct);
        for f=1:length(fs)
            sub=mystruct.(fs{f});
            if (~iscell(sub))
                sub={sub};
            end
            for ind=1:length(sub)
                if isstruct(sub{ind})
                    if isfield(values,fs{f})
                        check_field = values.(fs{f});
                        if isstruct(check_field) || iscell(check_field)
                            for structind=1:length(values.(fs{f}))
                                validatepaths(aap,mystruct.(fs{f})(1),values.(fs{f})(structind),[nme '.' fs{f}]);
                            end
                        else
                            validatepaths(aap,mystruct.(fs{f})(1),values.(fs{f}),[nme '.' fs{f}]);
                        end
                    end
                else
                    aas_log(aap,true,sprintf('Definition of field %s.%s in XML does not contain any descriptions (e.g., ui)\n',nme,fs{f}));
                end
            end
        end
    end
end

end
