%% Recursively copy across aap parameters from structure
function destaap=aas_copyparameters(srcaap,destaap,nme)

fn=fieldnames(srcaap);
for fnind=1:length(fn)
    if (isstruct(srcaap.(fn{fnind})))
        destaap.(fn{fnind})=aas_copyparameters(srcaap.(fn{fnind}),destaap.(fn{fnind}),[nme '.' fn{fnind}]);
    else
        if (~isfield(destaap,fn{fnind}))
            aas_log(srcaap,true,sprintf('Error when copying extra parameters, field %s is present in %s of extraparameters.aap  but not in normal aap structure',fn{fnind},nme));
        end;
        destaap.(fn{fnind})=srcaap.(fn{fnind});
    end;
end;
end
