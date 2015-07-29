function [aap, resp] = aamod_oneway_ANOVA(aap, task)
%
% -- conorwild -- 2014-05-02

resp='';

switch task
    case 'report'
        
    case 'doit'
        
        settings = aap.tasklist.currenttask.settings;
        UFp = 0.001;
        
        if ischar(settings.cells)
            settings.cells = str2num(settings.cells);
        end
        
        % New option to allow suffix to output file in extraparameters
        if (isfield(aap.tasklist.currenttask.extraparameters,'stats_suffix'))
            stats_suffix=aap.tasklist.currenttask.extraparameters.stats_suffix;
        else
            stats_suffix=[];
        end;
        
        % make analysis directory
        rfxrootdir = fullfile(aap.acq_details.root,[aap.directory_conventions.rfx stats_suffix]);
        if ~exist(rfxrootdir,'file'); mkdir(aap.acq_details.root,[aap.directory_conventions.rfx stats_suffix]);end
        cd(rfxrootdir);
        
        rfxdir = fullfile(rfxrootdir, settings.name);
        if exist(rfxdir)~=7; mkdir(rfxrootdir, settings.name);end
        cd(rfxdir);
        
        nsub = length(aap.acq_details.subjects);
        ncon = length(settings.cells);
        ncol = ncon; % Could be > if we use multiple BFs?
        nscan = sum(ncon .* nsub);
        conFiles = {};
        
        % Load SPMs and contrasts for each subject
        for subI = 1 : nsub
            spmFile = aas_getfiles_bystream(aap, subI, 'firstlevel_spm');
            SPM(subI) = load(spmFile);
            
            conFiles{subI} = aas_getfiles_bystream(aap, subI, 'firstlevel_cons');
            
        end
        

        
        allCellImgs = [];
        colNames = {};
        for conI = settings.cells
            
            % Make sure it's the same contrast across all subjects...
            conNames =  arrayfun(@(x) x.SPM.xCon(conI).name, SPM, 'UniformOutput', false);
            [conNames, subI] = unique(conNames);
            if length(conNames) > 1
                aas_log(aap, 1, sprintf('Contrast names are inconsistent across subjects! con_%05d returns: ''%s''', conI, strjoin(conNames, ', ')));
            end
            colNames = [colNames conNames];
            
            for sub = 1 : nsub
                origConFileName = SPM(sub).SPM.xCon(conI).Vcon.fname;
                matchingFile = regexp(cellstr(conFiles{sub}), origConFileName);
                matchingFile = ~cellfun(@(x) isempty(x), matchingFile);
                if ~any(matchingFile), error('Can''t find matching contrast file... more info added here eventuall'); end
                if sum(matchingFile)>1, error('More than one matching contrast image for this subject'); end
                allCellImgs = char(allCellImgs, conFiles{sub}(matchingFile,:));
            end

        end
        
        % Pop the 1st empty name.
        allCellImgs(1,:) = [];
        
        % Add in subject names (there are subject columns!)
        colNames = [colNames {aap.acq_details.subjects.mriname}];
        
        
        %-Assemble SPM structure
        %=======================================================================
        clear SPM;
        SPM.nscan = nscan;
        SPM.xY.P = allCellImgs;
        SPM.xY.VY = spm_vol(allCellImgs);
        
        % Build design matrix (X), Indices (Ind) and NONSPHERICITY (vi)...
        % Inelegant, but gets there...
        
        g = 1; ngrp = 1;
        X=[]; Ind=[]; vi={};
        nv=0; z=zeros(nscan,nscan);
        
        nr = ncol*nsub;
        id = [1:nsub]';
        
        X = [X; zeros(nr,ncol*(g-1)),...
            kron(eye(ncol), ones(nsub,1)),...
            zeros(nr,ncol*(ngrp-g))];
        
        % could add constants for group effects if wish
        
        % adding subject specific effects... not sure how this should
        % affect the sphericity correction
        % and this is probably wrong for between subjects ANOVA
        
        if settings.subjectfactor
            X = [X kron(ones(ncol,1),eye(nsub))];
            subjFactor = kron(ones(ncol,1),(1:nsub)');
        else
            subjFactor = ones(nscan,1);
        end
        
        SPM.factor.variance = settings.variance;
        SPM.factor.dept = settings.dependent;
        
        % The factor strucure Indices, for calculating non-sphericities
        factor1Index = kron([1:ncon]', ones(nsub,1));
        factor2Index = ones(nscan, 1);
        
        Ind = [ones(nscan, 1) factor1Index factor2Index  subjFactor ones(nscan,1)];
        
        SPM.xX = struct(...
            'X',X,...
            'iH',[1:(ncon)],'iC',zeros(1,0),'iB',zeros(1,0),'iG',zeros(1,0),...
            'name',{colNames},'I', Ind,...
            'sF',{{'repl'  'Contrast' 'dummy1' 'dummy2'  'grp'}});
        
        SPM.xVi.I = Ind;
        SPM = spm_get_vc(SPM);
        
        SPM.xC = [];
        
        SPM.xGX = struct(...
            'iGXcalc',1,    'sGXcalc','omit',                               'rg',[],...
            'iGMsca',9,     'sGMsca','<no grand Mean scaling>',...
            'GM',0,         'gSF', ones(nscan,1),...
            'iGC',  12,     'sGC', '(redundant: not doing AnCova)',        'gc',[],...
            'iGloNorm',9,   'sGloNorm','<no global normalisation>');
        
        %SPM.xVi = struct(...
        %        'iid',0,                'I',SPM.xX.I,                   'sF','SPM.xX.sF',...
        %        'var',[0 1 0 0],        'dep',[0 1 0 0],                'Vi',{vi} );
        
        % random modifications...
        % SPM.xVi = struct(...
        %     'I',SPM.xX.I, ...
        %     'Vi',{vi} );
        
        
        % Could add masking here...
        % set up masking here...
        
        Mdes    = struct(...
            'Analysis_threshold',   {'None (-Inf)'},...
            'Implicit_masking',     {'Yes: NaNs treated as missing'},...
            'Explicit_masking',     {'No'});
        
        SPM.xM  = struct(...
            'T',-Inf,'TH',ones(nscan,1)*-Inf,...
            'I',1,'VM',[],'xs',Mdes);
        
        Pdes = {{sprintf('%d condition, +0 covariate, +0 block, +0 nuisance',ncon); sprintf('%d total, having %d degrees of freedom',ncon,ncon); sprintf('leaving %d degrees of freedom from %d images',nscan-ncon,nscan)}};
        
        SPM.xsDes = struct(...
            'Design',               {'2-way ANOVA'},...
            'Global_calculation',   {'omit'},...
            'Grand_mean_scaling',   {'<no grand Mean scaling>'},...
            'Global_normalisation', {'<no global normalisation>'},...
            'Parameters',           Pdes);
        
        SPM.SPMid       = 'SPM8: (poorly) hacked together by Conor. ';
        
        save SPM SPM
        
        %save spmdes SPM
        
        % Estimate parameters
        %===========================================================================
        SPM = spm_spm(SPM);
        
        
        % Calculate the main effect contrast
        Ccon = ones(ncon,1)
        Dcon = orth(diff(eye(size(Ccon,1)))')';
        if settings.subjectfactor
            c = [Dcon zeros(size(Dcon,1), nsub)];
        else
            c = Dcon;
        end
        SPM.xCon = spm_FcUtil('Set', settings.name, 'F', 'c', c', SPM.xX.xKXs); 
        
        % Then compute the contrast to give average condition betas
        if settings.subjectfactor
            c = [eye(ncon) zeros(ncon,nsub)];
        else
            c = eye(ncon);
        end
        SPM.xCon(end+1) = spm_FcUtil('Set', 'Coondition Betas', 'F', 'c', c', SPM.xX.xKXs);   
        
        spm_contrasts(SPM);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
        
end

