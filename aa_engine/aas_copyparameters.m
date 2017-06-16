%% Recursively copy across aap parameters from structure
function destaap=aas_copyparameters(srcaap,destaap,nme)

fn=fieldnames(srcaap);
for fnind=1:length(fn)
    % index first struct element here to handle struct arrays
    if isstruct(srcaap(1).(fn{fnind}))
        destaap.(fn{fnind})=aas_copyparameters(srcaap.(fn{fnind}),destaap.(fn{fnind}),[nme '.' fn{fnind}]);
    else
        % index first struct element here to handle struct arrays
        if ~isfield(destaap(1),fn{fnind})
            aas_log(srcaap,true,sprintf('Error when copying extra parameters, field %s is present in %s of extraparameters.aap  but not in normal aap structure',fn{fnind},nme));
        end
        % struct array handling
        for n = 1:numel(srcaap)
            destaap(n).(fn{fnind})= srcaap(n).(fn{fnind});
        end
    end
end
