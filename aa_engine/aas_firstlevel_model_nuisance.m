function [movementRegs, compartmentRegs, physiologicalRegs, spikeRegs, GLMDNregs] = ...
    aas_firstlevel_model_nuisance(aap, subj, files)

settings = aap.tasklist.currenttask.settings;
streamIn = aap.tasklist.currenttask.inputstreams.stream;

% We can now have missing sessions per subject, so we're going to use only
% the sessions that are common to this subject and selected_sessions
[numSess, sessInds] = aas_getN_bydomain(aap, 'session', subj);
subjSessionI = intersect(sessInds, aap.acq_details.selected_sessions);
numSess = numel(subjSessionI);

for s = 1 : length(streamIn)
   if isstruct(streamIn{s})
       streamIn{s} = streamIn{s}.CONTENT;
   end
end

%% Movement regressors (extended!) [AVG]
movementRegs = [];

if isfield(aap.tasklist.currenttask.settings, 'includemovementpars') && ...
        aap.tasklist.currenttask.settings.includemovementpars == 1 && aas_stream_has_contents(aap,'realignment_parameter')
    
    [moves, mnames] = aas_movPars(aap,subj, aap.tasklist.currenttask.settings.moveMat);
    
    for sess = 1:numSess
        movementRegs(sess).names = mnames;
        movementRegs(sess).regs = moves{sess};
    end
end

%% Compartment regressors [AVG]
CRnames = {'GM', 'WM', 'CSF', 'OOH'};

compartmentRegs = [];

CRstreams = streamIn(strcmp('compSignal', streamIn));
if ~isempty(CRstreams) && aas_stream_has_contents(aap, 'compSignal');
    for sess=1:numSess
        % If we don't have compartment Signals, this should give up...
        compTC = [];
        load(aas_getfiles_bystream(aap,subj,subjSessionI(sess),CRstreams{:}));
        compartmentRegs(sess).names = CRnames(aap.tasklist.currenttask.settings.compRegs);
        compartmentRegs(sess).regs = compTC(:,aap.tasklist.currenttask.settings.compRegs);
    end
end

%% Physiological Regressors?
physiologicalRegs = [];

PRstreams = streamIn(strcmp('physreg', streamIn));
if ~isempty(PRstreams) && aas_stream_has_contents(aap, 'physreg')
    for sess=1:numSess
        PRfn = aas_getimages_bystream(aap,subj,subjSessionI(sess),PRstreams{:});
        
        % Contains physiological regressors
        R = []; names = [];
        load(PRfn);
        
        physiologicalRegs(sess).names = names;
        physiologicalRegs(sess).regs = R;
        
        if isempty(R)
            aas_log(aap,false, sprintf('Could not find Physiological Regressors for session %d\n', subjSessionI(sess)))
        end
    end
end

%% Spikes and moves, if these exist...
spikeRegs = [];

SPstreams = streamIn(strcmp('listspikes', streamIn));

if ~isempty(SPstreams) && aas_stream_has_contents(aap, 'listspikes') && settings.includespikes
	
    for sess = 1:numSess
		
        SPfn = aas_getimages_bystream(aap,subj,subjSessionI(sess),SPstreams{:});
        
        % Contains spike scan numbers
		
        TSspikes = []; Mspikes = [];
        load(SPfn);
        
		if (isempty(TSspikes))
			TSspikes = Mspikes;
		end

		if (~isempty(TSspikes))

			% Combine spikes and moves...
			
			regrscans = union(TSspikes(:,1), Mspikes(:,1));

			spikeRegs(sess).names = {};
			spikeRegs(sess).regs = zeros(size(files{sess},1), length(regrscans));

			for r = 1:length(regrscans)
				spikeRegs(sess).regs(regrscans(r),r) = 1;    % scan regrscan(r) is at one for scan r
				spikeRegs(sess).names{r} = sprintf('SpikeMov%d', r);
			end
			
		end
	
end
	


end

%% GLMdenoise regressors [CW]
GLMDNregs = [];
GLMDNstream = streamIn(strcmp('gd_results', streamIn));
if ~isempty(GLMDNstream) && aas_stream_has_contents(aap, 'gd_results') && settings.includeGLMDNregs
    GLMDNfile = aas_getfiles_bystream(aap, subj, GLMDNstream{:});
    load(GLMDNfile);
    for sess = 1:numSess
        GLMDNregs(sess).names = arrayfun(@(x) sprintf('GLMDN_%02d', x), [1:gd_results.pcnum], 'UniformOutput', false);
        GLMDNregs(sess).regs = gd_results.pcregressors{sess}(:, 1:gd_results.pcnum);
    end
end

