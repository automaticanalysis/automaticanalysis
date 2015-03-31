% Run a group ICA analysis with 'gift'

function [aap, resp] = aamod_secondlevel_GIFT(aap, task)

resp='';

switch task
    case 'report'
        
    case 'doit'
        settings = aap.tasklist.currenttask.settings;
        localroot = aas_getstudypath(aap);
        
        %% Paths
        load GIFT_Defaults
        sesInfo.userInput.prefix = 'gift';
		sesInfo.userInput.pwd = fullfile(localroot,'gift');
        if ~exist(sesInfo.userInput.pwd,'dir'), mkdir(sesInfo.userInput.pwd); end
		sesInfo.userInput.param_file = fullfile(sesInfo.userInput.pwd,...
            [sesInfo.userInput.prefix '_ica_parameters.mat']);
        
        %% Data
        subjects = aap.acq_details.subjects;
        if ~isempty(settings.selected_subjects), subjects = subjects(settings.selected_subjects); end
        sessions = aap.acq_details.sessions;
        if ~isempty(settings.selected_session), sessions = sessions(settings.selected_session); end
        for subj = 1:numel(subjects)
            firstlevel_brainmask(subj,:) = aas_getfiles_bystream(aap,'subject',subj,'firstlevel_brainmask');
            for sess = 1:numel(sessions)
                fepi = aas_getfiles_bystream(aap,'session',[subj sess],'epi');
                if aap.options.NIFTI4D
                    fepi = spm_select('ExtFPList',fileparts(fepi),basename(fepi));
                end
                ind = (subj-1)*numel(sessions)+sess;
                sesInfo.userInput.files(ind) = struct('name', fepi);
                sesInfo.userInput.diffTimePoints(ind) = size(fepi,1);
            end
        end
        sesInfo.userInput.numOfSub = subj;
        sesInfo.userInput.numOfSess = sess;
        sesInfo.userInput.numOfGroups1 = sesInfo.userInput.numOfSub*sesInfo.userInput.numOfSess;

        mask = spm_read_vols(spm_vol(firstlevel_brainmask));
        mask = min(mask,[],4);
        V = spm_vol(firstlevel_brainmask(1,:));        
        sesInfo.userInput.maskFile = fullfile(sesInfo.userInput.pwd,'Mask.nii');
        nifti_write(sesInfo.userInput.maskFile,mask,'GIFT Mask',V);
        V = spm_vol(sesInfo.userInput.files(1).name(1,:));
        [junk,VOX] = spm_get_bbox(V); VOX = abs(VOX);
        sesInfo.userInput.HInfo = struct(...
            'V',V,...
            'DIM',V.dim,...
            'VOX',VOX...
            );
        sesInfo.userInput.mask_ind = find(mask);
        
        %% Estimate dimensions
        if settings.dimest
            mask = zeros(sesInfo.userInput.HInfo.DIM(1:3));
            mask(sesInfo.userInput.mask_ind) = 1;
            mask = (mask == 1);
            sessind = sesInfo.userInput.files;
            if ~isempty(settings.subject_to_dimest)
                sessind = settings.subject_to_dimest;
            end
            for sess = 1:numel(sessind)
                sesInfo.userInput.estimated_comps(sess) = icatb_estimate_dimension(sesInfo.userInput.files(sess).name,mask,icatb_checkPrecision('single', 0));
            end
        end
        
         %% ICA parameters
        if isfield(sesInfo.userInput,'estimated_comps')
            if ~ischar(settings.numICs)
                aas_log(aap,true,'Dimensionality is estimated, therfore numICs should specify method (e.g. mean, min, max)');
            end
            sesInfo.userInput.numComp = round(eval([lower(settings.numICs), '(sesInfo.userInput.estimated_comps);']));
            sesInfo.userInput.numOfPC1 = round(1.5*sesInfo.userInput.numComp);
            sesInfo.userInput.numOfPC2 = sesInfo.userInput.numComp;
        else
            if isempty(settings.numICs) || ~isnumeric(settings.numICs)
                aas_log(aap,true,'Number of components must be specified!');
            else
                sesInfo.userInput.numComp = settings.numICs;
            end
            if ~isempty(settings.PCA1)
                sesInfo.userInput.numOfPC1 = settings.PCA1;
            else
                sesInfo.userInput.numOfPC1 = round(1.5*sesInfo.userInput.numComp);
            end
            if ~isempty(settings.PCA2)
                sesInfo.userInput.numOfPC2 = settings.PCA2;
            else
                sesInfo.userInput.numOfPC2 = sesInfo.userInput.numComp;
            end
        end

        sesInfo.userInput.algorithm = settings.algorithm;
        sesInfo.userInput.preproc_type = settings.preproc;
        sesInfo.userInput.ICA_Options = icatb_icaOptions([sesInfo.userInput.numComp, length(sesInfo.userInput.mask_ind)],sesInfo.userInput.algorithm,'off');
        
        sesInfo.userInput.scaleType = settings.scaleType;
        
        switch settings.analysis
            case 2
                sesInfo.userInput.icasso_opts = settings.icasso;
            otherwise
        end
        
        %% Save parameters ans run
        save(sesInfo.userInput.param_file,'sesInfo');
        icatb_runAnalysis(sesInfo, 1)
        
        %% Outputs
        aap = aas_desc_outputs(aap, 'study', [], 'gift_parameters', fullfile(sesInfo.userInput.param_file));
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end

end