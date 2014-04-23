function [movementRegs, compartmentRegs, physiologicalRegs, spikeRegs, GLMDNregs] = ...
    aas_firstlevel_model_nuisance(aap, subj, files)

streamIn = aap.tasklist.currenttask.inputstreams.stream;

%% Movement regressors (extended!) [AVG]
movementRegs = [];

if isfield(aap.tasklist.currenttask.settings, 'includemovementpars') && ...
        aap.tasklist.currenttask.settings.includemovementpars == 1 && aas_stream_has_contents(aap,'realignment_parameter')
    
    [moves, mnames] = aas_movPars(aap,subj, aap.tasklist.currenttask.settings.moveMat);
    
    for sess = aap.acq_details.selected_sessions
        movementRegs(sess).names = mnames;
        movementRegs(sess).regs = moves{sess};
    end
end

%% Compartment regressors [AVG]
CRnames = {'GM', 'WM', 'CSF', 'OOH'};

compartmentRegs = [];

CRstreams = streamIn(strcmp('compSignal', streamIn));
if ~isempty(CRstreams)
    for sess=aap.acq_details.selected_sessions
        % If we don't have compartment Signals, this should give up...
        compTC = [];
        load(aas_getfiles_bystream(aap,subj,sess,CRstreams{:}));
        compartmentRegs(sess).names = CRnames(aap.tasklist.currenttask.settings.compRegs);
        compartmentRegs(sess).regs = compTC(:,aap.tasklist.currenttask.settings.compRegs);
    end
end

%% Physiological Regressors?
physiologicalRegs = [];

PRstreams = streamIn(strcmp('physreg', streamIn));
if ~isempty(PRstreams)
    for sess=aap.acq_details.selected_sessions
        PRfn = aas_getimages_bystream(aap,subj,sess,PRstreams{:});
        
        % Contains physiological regressors
        R = []; names = [];
        load(PRfn);
        
        physiologicalRegs(sess).names = names;
        physiologicalRegs(sess).regs = R;
        
        if isempty(R)
            aas_log(aap,false, sprintf('Could not find Physiological Regressors for session %d\n', sess))
        end
    end
end

%% Spikes and moves, if these exist...
spikeRegs = [];

SPstreams = streamIn(strcmp('listspikes', streamIn));

if ~isempty(SPstreams)
    for sess=aap.acq_details.selected_sessions
        SPfn = aas_getimages_bystream(aap,subj,sess,SPstreams{:});
        
        % Contains spike scan numbers
        TSspikes = []; Mspikes = [];
        load(SPfn);
        
        % Combine spikes and moves...
        regrscans = union(TSspikes(:,1), Mspikes(:,1));
        
        spikeRegs(sess).names = {};
        spikeRegs(sess).regs = zeros(size(files{sess},1), length(regrscans));
        
        for r = 1:length(regrscans),
            spikeRegs(sess).regs(regrscans(r),r) = 1;    % scan regrscan(r) is at one for scan r
            spikeRegs(sess).names{r} = sprintf('SpikeMov%d', r);
        end
    end
end

%% GLMdenoise regressors [CW]
GLMDNregs = [];
GLMDNstream = streamIn(strcmp('gd_results', streamIn));
if ~isempty(GLMDNstream)
    GLMDNfile = aas_getfiles_bystream(aap, subj, GLMDNstream{:});
    load(GLMDNfile);
    for sess = aap.acq_details.selected_sessions
        GLMDNregs(sess).names = arrayfun(@(x) sprintf('GLMDN_%02d', x), [1:gd_results.pcnum], 'UniformOutput', false);
        GLMDNregs(sess).regs = gd_results.pcregressors{sess}(:, 1:gd_results.pcnum);
    end
end