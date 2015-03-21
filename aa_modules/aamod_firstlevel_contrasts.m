% AA module
% First-level model Adam Hampshire MRC CBU Cambridge Feb 2006
% Modified for aa by Rhodri Cusack Mar 2006-2011
% Additions by Rik Henson Mar 2011
% Tibor Auer MRC CBU Cambridge 2012-2013

function [aap,resp]=aamod_firstlevel_contrasts(aap,task,subj)

resp='';


switch task
    case 'domain'
        resp='subject'; % this module needs to be run once per subject
        
    case 'description'
        resp='SPM5 contrasts';
        
    case 'summary'
        subjpath=aas_getsubjpath(subj);
        resp=sprintf('Contrasts %s\n',subjpath);
        
    case 'report' % [TA]
        if ~exist(fullfile(aas_getsubjpath(aap,subj),'diagnostic_aamod_firstlevel_contrast.jpg'),'file')
            efficiency(aap,subj);
        end
        fdiag = dir(fullfile(aas_getsubjpath(aap,subj),'diagnostic_*.jpg'));
        for d = 1:numel(fdiag)
            aap = aas_report_add(aap,subj,'<table><tr><td>');
            aap=aas_report_addimage(aap,subj,fullfile(aas_getsubjpath(aap,subj),fdiag(d).name));
            aap = aas_report_add(aap,subj,'</td></tr></table>');
        end
        
    case 'doit'
        mriname = aas_prepare_diagnostic(aap,subj);
        
        cwd=pwd;
        % get the subdirectories in the main directory
        subj_dir  =  aas_getsubjpath(aap,subj);
        
        % Maintained for backwards compatibility- better now now put
        % module-specific value in
        % aap.directory_conventions.stats_singlesubj
        if (isfield(aap.tasklist.currenttask.extraparameters,'stats_suffix'))
            stats_suffix=aap.tasklist.currenttask.extraparameters.stats_suffix;
        else
            stats_suffix=[];
        end
        anadir  =  fullfile(subj_dir,[aap.directory_conventions.stats_singlesubj stats_suffix]);
        
        % Now set up contrasts...
        SPM=load(aas_getfiles_bystream(aap,subj,'firstlevel_spm'));
        SPM=SPM.SPM;
        SPM.swd=anadir;
        
        settings = aap.tasklist.currenttask.settings;
        
        if settings.useaasessions
            [nsess, sessInds] = aas_getN_bydomain(aap, 'session', subj);
            subjSessionI = intersect(sessInds, aap.acq_details.selected_sessions);
            sessnames = {aap.acq_details.sessions(subjSessionI).name};
            selected_sessions = subjSessionI;
            nsess = length(selected_sessions);
            nsess_all = length(sessnames);
        else
            % just get all sessions based on SPM file
            sessnames = {};
            nsess = length(SPM.Sess);
            selected_sessions = 1:nsess;
            nsess_all = nsess;
        end        
        
        % Load up contrasts from task settings
        [fle, subjname, ext] = fileparts(subj_dir);
        contrasts_set = find(strcmp({settings.contrasts.subject}, [subjname ext]));
        if (isempty(contrasts_set))
            % Try for wildcard
            contrasts_set=find(strcmp({settings.contrasts.subject},'*'));
            if (isempty(contrasts_set))
                aas_log(aap,true,'Can''t find declaration of what contrasts to use  -  insert this in a local copy of aamod_firstlevel_contrasts.xml or put into user script');
            end
        end
        
        contrasts=settings.contrasts(contrasts_set);
        % add contrasts for each task regressor v baseline?
        if settings.eachagainstbaseline
            basev = zeros(1,length(SPM.Sess(1).col));
            for conind = 1:length(basev)
                newv = basev;
                newv(conind) = 1;
                contrasts.con(end+1)= struct('name',sprintf(...
                    '%s-o-baseline',SPM.xX.name{SPM.xX.iC(conind)}),...
                    'format','sameforallsessions',...
                    'vector',newv,...
                    'session',[],...
                    'contype','T');
            end
        end
        
        % logical vector for run-specific contrasts
        % First the general case - across all runs
        nregr = length(SPM.xX.name);
        nruns = length(SPM.Sess);
        runI = true(1,nregr);
        noregr = zeros(1,nregr);
        
        % Do we also want run-specific contrasts?
        if settings.oneconperrun && (nruns > 1)
            for r = 1:nruns
                % All zeros
                runI(r+1,:) = noregr;
                % Then fill in run regs
                runI(r+1,SPM.Sess(r).col)=1;
            end
        end
        
        ccount = 0;
        convec_names  = SPM.xX.name(SPM.xX.iC);
        
        % Separately for each run config (only one if ~oneconperrun)
        for r = 1:size(runI,1)
            % For each unique contrast
            for conind=1:length(contrasts.con)
                % For each individual contrast (different if oneconperrun)
                ccount = ccount + 1;
                % Get or make automatic contrast name
                if (isfield(contrasts.con(conind),'name') && ~isempty(contrasts.con(conind).name))
                    finalname = contrasts.con(conind).name;
                else
                    finalname = sprintf('Con%d',conind);
                end
                % may have to change name to reflect run
                if r == 1
                    connames{ccount} = finalname;
                else
                    connames{ccount} = sprintf('%s-run%02d',finalname,r-1);
                end
                % support eval'ed strings to define contrasts (e.g. ones, eye)
                if ischar(contrasts.con(conind).vector)
                    contrasts.con(conind).vector = eval(...
                        contrasts.con(conind).vector);
                end
                % Make contract vector
                switch(contrasts.con(conind).format)
                    
                    case {'singlesession','sameforallsessions'}
                        if (strcmp(contrasts.con(conind).format,'singlesession'))
                            sessforcon=[strcmp(sessnames,contrasts.con(conind).session)];
                        else
                            % [AVG] To make the selected sessions work...
                            sessforcon = zeros(1,nsess_all);
                            for sess = selected_sessions
                                sessforcon(sess) = 1;
                            end
                        end
                        convec=[];
                        sessnuminspm=1;
                        for sess=selected_sessions
                            numcolsinthissess = length(SPM.Sess(sessnuminspm).col);
                            if (sessforcon(sess))
                                if (size(contrasts.con(conind).vector,2) > numcolsinthissess)
                                    aas_log(aap,true,sprintf('Number of columns in contrast matrix for session %d is more than number of columns in model for this session - wanted %d columns, got ',sess,numcolsinthissess)); disp(contrasts.con(conind).vector);
                                elseif (size(contrasts.con(conind).vector,2) < numcolsinthissess)
                                    convec = [convec contrasts.con(conind).vector zeros(size(contrasts.con(conind).vector,1),numcolsinthissess-size(contrasts.con(conind).vector,2))];
                                    %aas_log(aap,false,sprintf('Warning: Number of columns in contrast matrix for session %d is less than number of columns in model for this session - wanted %d columns, so padding to ',sess,numcolsinthissess)); disp(convec);
                                else
                                    convec = [convec contrasts.con(conind).vector];
                                end
                            else
                                convec = [convec zeros(size(contrasts.con(conind).vector,1), numcolsinthissess)];
                            end
                            sessnuminspm=sessnuminspm+1;
                        end
                        
                        % If subjects have different # of sessions, then
                        % they will be weighted differently in 2nd level
                        % model. So, normalize the contrast by the number
                        % of sesisons that contribute to it [CW]
                        convec = convec ./ nnz(sessforcon); 
                        
                    case 'uniquebysession'
                        totnumcolsbarconstants = size(SPM.xX.X,2) - nsess;
                        
                        if (size(contrasts.con(conind).vector,2) > totnumcolsbarconstants)
                            aas_log(aap,true,sprintf('Number of columns in contrast matrix for session %d is more than number of columns in model (bar constants) - wanted %d columns, got ',totnumcolsbarconstants)); disp(contrasts.con(conind).vector);
                        elseif (size(contrasts.con(conind).vector,2) < totnumcolsbarconstants)
                            convec = contrasts.con(conind).vector;
                            if settings.automatic_movesandmeans
                                % [AVG] *better* way of specifying the correct columns...
                                convec_out = zeros(1,totnumcolsbarconstants);
                                convec_out(SPM.xX.iC) = convec;
                                                                
                                convec = convec_out;
                            end
                        else
                            convec=contrasts.con(conind).vector;
                        end
                    otherwise
                        aas_log(aap,true,sprintf('Unknown format %s specified for contrast %d',contrasts.con(conind).format,ccount));
                end
                cons{ccount} = [convec zeros(size(convec,1),nsess)];  % Add final constant terms
                
                % Check not empty
                if (~any(cons{ccount}(:)))
                    aas_log(aap,true,sprintf('Contrast %d has no non-zero values, not permitted.',contrasts_set(ccount)));
                end
                
                % Allow F tests
                if (isfield(contrasts.con(conind),'contype') && isempty(contrasts.con(conind).contype))
                    contype{ccount}='T';
                else
                    contype{ccount}=contrasts.con(conind).type;
                end
                
                % Zero out run-irrelevant entries
                % support for multi-row F contrasts
                nrows = size(cons{ccount},1);
                inds = repmat(runI(r,:)~=1,[nrows 1]);
                cons{ccount}(inds) = 0;
                
                % DIAGNOSTIC
                fprintf('\n%s\n', contrasts.con(conind).name)
                for conind  =  1:max(size(convec_names))
                    fprintf('\t%s: %d\n', convec_names{conind}, convec(SPM.xX.iC(conind)))
                end
            end
        end
        
        % Make the con images
        SPM.xCon =[];
        for conind = 1:length(cons)
            % skip empty regressors
            if all(cons{conind}(:) == 0)
                continue
            end
            if length(SPM.xCon)==0
                SPM.xCon = spm_FcUtil('Set', connames{conind}, contype{conind},'c', cons{conind}', SPM.xX.xKXs);
            else
                SPM.xCon(end+1) = spm_FcUtil('Set', connames{conind}, contype{conind},'c', cons{conind}', SPM.xX.xKXs);
            end
        end
        SPM = spm_contrasts(SPM);
        
        % Efficiency based on Rik Henson's script [TA]
        if settings.estimateefficiency
            efficiency(aap, subj, SPM); 
        end
        
        % Describe outputs
        %  updated spm
        aap = aas_desc_outputs(aap,subj,'firstlevel_spm',fullfile(anadir,'SPM.mat'));
        
        %  firstlevel_betas (includes related statistical files)
        filters={'con','spmT','spmF'};
        
        for filterind=1:length(filters)
            allbetas=dir(fullfile(anadir,[filters{filterind} '_*.nii']));
            betafns=[];
            for betaind=1:length(allbetas);
                betafns=strvcat(betafns,fullfile(anadir,allbetas(betaind).name));
            end
            allbetas=dir(fullfile(anadir,[filters{filterind} '_*.img']));
            for betaind=1:length(allbetas);
                betafns=strvcat(betafns,fullfile(anadir,allbetas(betaind).name));
            end
            allbetas=dir(fullfile(anadir,[filters{filterind} '_*.hdr']));
            for betaind=1:length(allbetas);
                betafns=strvcat(betafns,fullfile(anadir,allbetas(betaind).name));
            end
            aap=aas_desc_outputs(aap,subj,['firstlevel_' lower(filters{filterind}) 's'],betafns);
        end
        cd (cwd);
        
        %% DIAGNOSTICS (check distribution of T-values in contrasts)
%         D = dir(fullfile(anadir, 'spmT_*.img'));
%         for d = 1:length(D)
%             h = img2hist(fullfile(anadir, D(d).name), [], contrasts.con(d).name);
%             saveas(h, fullfile(aap.acq_details.root, 'diagnostics', ...
%                 [mfilename '__' mriname '_' contrasts.con(d).name '.fig']), 'fig');
%             try close(h); catch; end
%         end
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
        
end
end

function h = efficiency(aap,subj,SPM)
% Based on Rik Henson's script

% Note this calculation of efficiency takes the 'filtered and whitened'
% design matrix (X) as it is in SPM.

% SPM gets passed in now, because this will load the input stream SPM file,
% which doesn't (always) have the contrasts.

if nargin < 3 % SPM is not passed (e.g. reporting)
    load(aas_getfiles_bystream(aap,subj,'firstlevel_spm'));
end
X = SPM.xX.xKXs.X;
iXX=inv(X'*X);

[junk, nameCols] = strtok(SPM.xX.name(SPM.xX.iC),' ');
nameCols = strtok(nameCols,'*');
nameCons = {SPM.xCon.name}';
cons = {SPM.xCon.c};
effic = nan(numel(cons), 1);
columnsCon = nan(numel(cons), length(nameCols));

for conind = 1:numel(cons)    
    fullCon = cons{conind}';
    selCon = fullCon(SPM.xX.iC);    
    columnsCon(conind, :) = selCon;    
    
    % Normalize Contrast
    fullCon = fullCon / max(sum(fullCon(fullCon>0)), sum(fullCon(fullCon < 0)));
    % Calculate efficiency
    effic(conind) = trace(fullCon*iXX*fullCon')^-1;    
end

% get Text size
h = figure; set(h, 'Position', [1 1 1280 720]); % HD
ht = text(1,1,nameCons,'FontSize',12,'FontWeight','Bold','interpreter','none');
set(ht,'Unit','normalized');
tSize = get(ht,'Extent'); tWidth = tSize(3);
close(h);

h = figure; set(h, 'Position', [1 1 1280 720]); % HD
subplot('Position', [tWidth 0.1 0.6-tWidth 0.9/20*numel(cons)]); % assume not more then 20 contrast
imagesc(columnsCon)
set(gca, 'YTick', 1:numel(cons), 'YTickLabel',nameCons,  ...
    'Xtick', 1:length(nameCols), 'XTickLabel',nameCols)
set(gca, 'XAxisLocation','top');
xlab = rotateticklabel(gca,90);
set(gca,'FontSize',12,'FontWeight','Bold');
set(xlab,'FontSize',12,'FontWeight','Bold');

subplot('Position', [0.6 0.1 0.35 0.9/20*numel(cons)]);
set(gca, 'YTick', 1:numel(cons), 'YTickLabel','');
set(gca,'FontSize',12,'FontWeight','Bold');
hold on;
cmap = colorcube(numel(cons));

for conind = 1:numel(cons)        
    barh(numel(cons) - conind + 1, log(effic(conind)), 'FaceColor',cmap(conind,:));
end

ylim([0.5 numel(cons)+0.5])
xlabel('Log Efficiency')
efficiencyVals = floor(log(min(effic))):0.1:ceil(log(max(effic)));
set(gca, 'Xtick', efficiencyVals, 'XtickLabel', sprintf('%1.1f|',exp(efficiencyVals)))

fname = fullfile(aas_getsubjpath(aap,subj),'diagnostic_aamod_firstlevel_contrast.jpg');
set(h,'Renderer','zbuffer');
print(h,'-djpeg','-r150',fname);
close(h);

% cut image assuming (1,1) to be empty
cdata = imread(fname);
img = sum(cdata,3); rows = sum(img,2); columns = sum(img,1);
row1 = find(rows ~= rows(1),1,'first');
row2 = find(rows ~= rows(1),1,'last');
column1 = find(columns ~= columns(1),1,'first');
column2 = find(columns ~= columns(1),1,'last');
cdata = cdata(row1:row2,column1:column2,:);
imwrite(cdata,fname);
end
