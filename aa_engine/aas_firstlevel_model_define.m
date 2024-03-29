function [SPM, cols_interest, cols_nuisance, currcol] = aas_firstlevel_model_define(aap, sess, sessnuminspm, SPM, model, modelC, ...
    cols_interest, cols_nuisance, currcol, ...
    movementRegs, compartmentRegs, physiologicalRegs, spikeRegs, GLMDNregs)


%% Get Information about Basis functions 
%  could have more than one column:
%  e.g., time derivatives or a FLOBS-type basis set

xBF = SPM.xBF;

% Use an SPM BF if one is isn't specified
if ~isfield(xBF, 'bf') || isempty(xBF.bf)
    xBF = [];
    xBF.dt = SPM.xY.RT;
    xBF.name = SPM.xBF.name;
    xBF.length = SPM.xBF.length;
    xBF.order = SPM.xBF.order;
    xBF = spm_get_bf(xBF);
end

% How many regressors make up this BF?
numBFregs = size(xBF.bf, 2);

% Create empty session, in case we don't actually find a model or anything
U = struct( 'ons', {},...
            'dur', {},...
            'name', {{}},...
            'P', {},...
            'orth', [] );

SPM.Sess(sessnuminspm).U = U;

%% Define model{sess} events
if ~isempty(model{sess})
    for c = 1:length(model{sess}.event);
        if (isempty(model{sess}.event(c).parametric))
            parametric=struct('name','none');
            parLen = 0;
        else
            parametric=model{sess}.event(c).parametric;
            parLen = length(parametric);
            
            % scale temporal modulation
            tmods_to_scale = find(strcmp({parametric.name},'time_tosc'));
            if ~isempty(tmods_to_scale)
                switch SPM.xBF.UNITS
                    case 'secs'
                        sf    = 1 / 60;
                    case 'scans'
                        sf    = SPM.xY.RT / 60;
                    otherwise
                        aas_log(aap,true,'ERROR: Unknown unit "%s".',SPM.xBF.UNITS);
                end            
                for p = tmods_to_scale
                    parametric(p).P = parametric(p).P * sf;
                end
            end
        end
        
        SPM.Sess(sessnuminspm).U(c) = struct(...
            'ons', model{sess}.event(c).ons,...
            'dur', model{sess}.event(c).dur,...
            'name', {{model{sess}.event(c).name}},...
            'P',parametric,...
            'orth', 1 );
        
        % Orthogonalise
        if strcmp(spm('Ver'),'SPM12') &&...
                isfield(aap.tasklist.currenttask.settings,'orthogonalisation') &&...
                ~aap.tasklist.currenttask.settings.orthogonalisation
            SPM.Sess(sessnuminspm).U(c).orth = 0;
        end
        
        cols_interest=[cols_interest [1:(1+parLen)*numBFregs]+currcol-1];
        currcol=currcol+(1+parLen)*numBFregs;
    end
end 

%% Nuisance and covariate regressors
SPM.Sess(sessnuminspm).C.C = [];
SPM.Sess(sessnuminspm).C.name = {};

%% Define model covariates
if ~isempty(modelC{sess})
    %% Set up the convolution vector...
    % xBF.dt      - time bin length {seconds}
    % xBF.name    - description of basis functions specified
    % xBF.length  - window length (seconds)
    % xBF.order   - order

    for c = 1:length(modelC{sess}.covariate);
        covVect = modelC{sess}.covariate(c).vector;
        
        % Do we convolve with HRF?
        if modelC{sess}.covariate(c).HRF > 0
            U =[];
            U.u = covVect(:);
            U.name = {modelC{sess}.covariate(c).name};
            covVect = spm_Volterra(U, xBF.bf);
        end
        
        SPM.Sess(sessnuminspm).C.C    = [SPM.Sess(sessnuminspm).C.C ...
            covVect];     % covariate
        SPM.Sess(sessnuminspm).C.name = [SPM.Sess(sessnuminspm).C.name ...
            modelC{sess}.covariate(c).name];
        
        % Is the covariate of interest or nuisance
        if modelC{sess}.covariate(c).interest > 0
            cols_interest=[cols_interest currcol];
        else
            cols_nuisance=[cols_nuisance currcol];
        end
        currcol=currcol + 1;
    end
end

%% Movement regressors
if ~isempty(movementRegs)
    SPM.Sess(sessnuminspm).C.C    = [SPM.Sess(sessnuminspm).C.C ...
        movementRegs(sess).regs];
    SPM.Sess(sessnuminspm).C.name = [SPM.Sess(sessnuminspm).C.name ...
        movementRegs(sess).names];
    
    cols_nuisance = [cols_nuisance currcol:(currcol + ...
        length(movementRegs(sess).names)) - 1];
    currcol = currcol ...
        + length(movementRegs(sess).names);
end

%% Compartment regressors [AVG]
if ~isempty(compartmentRegs)
    SPM.Sess(sessnuminspm).C.C    = [SPM.Sess(sessnuminspm).C.C ...
        compartmentRegs(sess).regs];
    SPM.Sess(sessnuminspm).C.name = [SPM.Sess(sessnuminspm).C.name ...
        compartmentRegs(sess).names];
    
    cols_nuisance = [cols_nuisance currcol:(currcol + ...
        + length(aap.tasklist.currenttask.settings.compRegs) - 1)];
    currcol = currcol ...
        + length(aap.tasklist.currenttask.settings.compRegs);
end

%% Physiological Regressors?
if ~isempty(physiologicalRegs)
    SPM.Sess(sessnuminspm).C.C = [SPM.Sess(sessnuminspm).C.C ...
        physiologicalRegs(sess).regs];
    SPM.Sess(sessnuminspm).C.name = [SPM.Sess(sessnuminspm).C.name ...
        physiologicalRegs(sess).names];
    
    cols_nuisance=[cols_nuisance currcol:(currcol+length(physiologicalRegs(sess).names) - 1)];
    currcol = currcol + length(physiologicalRegs(sess).names);
end

%% Spikes and moves, if these exist...

% it's possible for spikeRegs to be defined for some sessions and not
% others (i.e., spikeRegs won't be empty but spikeRegs(sess) might be)
% so we have to do some additional checking here...

if (~isempty(spikeRegs) && le(sess,numel(spikeRegs)) && ~isempty(spikeRegs(sess).regs))
    SPM.Sess(sessnuminspm).C.C = [SPM.Sess(sessnuminspm).C.C spikeRegs(sess).regs];
    SPM.Sess(sessnuminspm).C.name = [SPM.Sess(sessnuminspm).C.name spikeRegs(sess).names];
    cols_nuisance=[cols_nuisance currcol:(currcol+length(spikeRegs(sess).names) - 1)];
    currcol = currcol + length(spikeRegs(sess).names);
end

%% GLMdenoise regressors, if they exist... [CW]
if ~isempty(GLMDNregs)
    SPM.Sess(sessnuminspm).C.C = [SPM.Sess(sessnuminspm).C.C ...
        GLMDNregs(sess).regs];
    SPM.Sess(sessnuminspm).C.name = [SPM.Sess(sessnuminspm).C.name ...
        GLMDNregs(sess).names];
    
    cols_nuisance=[cols_nuisance currcol:(currcol+length(GLMDNregs(sess).names) - 1)];
    currcol = currcol + length(GLMDNregs(sess).names);
end

end

