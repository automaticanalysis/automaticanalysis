% AA module
% First-level model Adam Hampshire MRC CBU Cambridge Feb 2006
% Modified for aa by Rhodri Cusack Mar 2006-2011
% Additions by Rik Henson Mar 2011
% Tibor Auer MRC CBU Cambridge 2012-2013
%
% CHANGE HISTORY
%
% 03/2023 [MSJ] fix eachagainstbaseline; remove 20 contrast restriction
%               on diagnostic plot; save a plaintext list of contrast names
% 09/2021 [MSJ] added string contrast option for "uniquebysession";

function [aap,resp]=aamod_firstlevel_contrasts(aap,task,subj)

resp='';


switch task
    case 'report' % [TA]
        if ~exist(fullfile(aas_getsubjpath(aap,subj),['diagnostic_' mfilename '.jpg']),'file')
            diag(aap,subj);
        end
        fdiag = dir(fullfile(aas_getsubjpath(aap,subj),'diagnostic_*.jpg'));
        for d = 1:numel(fdiag)
            aap = aas_report_add(aap,subj,'<table><tr><td>');
            aap=aas_report_addimage(aap,subj,fullfile(aas_getsubjpath(aap,subj),fdiag(d).name));
            aap = aas_report_add(aap,subj,'</td></tr></table>');
        end
        
    case 'doit'
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
            [~, sessInds] = aas_getN_bydomain(aap, 'session', subj);
            subjSessionI = intersect(sessInds, aap.acq_details.selected_sessions);
            sessnames = {aap.acq_details.sessions(:).name};
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
        contrasts_set = find(strcmp({settings.contrasts.subject}, basename(subj_dir)));
        if (isempty(contrasts_set))
            % Try for wildcard
            contrasts_set=find(strcmp({settings.contrasts.subject},'*'));
        end
        
        contrasts=settings.contrasts(contrasts_set);
        
        % add contrasts for each task regressor v baseline?
        if settings.eachagainstbaseline
            if isempty(contrasts), contrasts(1).subject = aas_getsubjname(aap,subj); end
            % contrast vector SHOULDNT include entries for the model
            % constants (causes crash later) - so just subtract off number
            % of sessions and we should be G2G (SPM.xX.iC should never
            % reference these slots) foobar
            basev = zeros(1,numel(SPM.xX.name)-nsess);
            for conind = 1:numel(SPM.xX.iC)
                newv = basev;
                newv(SPM.xX.iC(conind)) = 1; 
                contrasts.con(end+1)= struct(...
                    'format','uniquebysession',...
                    'vector',newv,...
                    'session',[],...
                    'type','T',...
                    'name',sprintf('%s-o-baseline',SPM.xX.name{SPM.xX.iC(conind)})...
                    );
            end
        end

        
        if isempty(contrasts)
            aas_log(aap,true,'ERROR: Can''t find declaration of what contrasts to use  -  check user master script!');
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
                if ischar(contrasts.con(conind).vector) &&...
                        ~any(strcmp({'+' '-'},contrasts.con(conind).vector(1))) % not defined with EV names
                    contrasts.con(conind).vector = str2num(...
                        contrasts.con(conind).vector);
                end
                % Make contrast vector
                switch(contrasts.con(conind).format)
                    
                    case {'singlesession','sessions','sameforallsessions'}
                        sessforcon = zeros(1,nsess_all);
                        if strcmp(contrasts.con(conind).format,'sameforallsessions')
                            % [AVG] To make the selected sessions work...
                            sessforcon(selected_sessions) = 1;
                        else
                            sessions = contrasts.con(conind).session;
                            if isstruct(sessions)
                                sessweights = sessions.weights;
                                sessions = sessions.names;                                
                            else
                                sessweights = ones(1,numel(sessions));
                            end
                            [~,indsess] = intersect(sessnames,sessions);
                            sessforcon(indsess) = sessweights;
                        end
                        
                        convec=[];
                        sessnuminspm=0;
                        for sess=selected_sessions
                            sessnuminspm=sessnuminspm+1;
                            numcolsinthissess = length(SPM.Sess(sessnuminspm).col);
                            if (sessforcon(sess))
                                if isnumeric(contrasts.con(conind).vector) % vcontrast vector
                                    if (size(contrasts.con(conind).vector,2) > numcolsinthissess)
                                        aas_log(aap,true,sprintf('ERROR: Number of columns in contrast matrix for session %d is more than number of columns in model for this session - wanted %d columns, got %d ',sess,numcolsinthissess,size(contrasts.con(conind).vector,2))); 
                                    elseif (size(contrasts.con(conind).vector,2) < numcolsinthissess) % padding if shorter
                                        convec = [convec sessforcon(sess)*contrasts.con(conind).vector zeros(size(contrasts.con(conind).vector,1),numcolsinthissess-size(contrasts.con(conind).vector,2))];
                                    else
                                        convec = [convec sessforcon(sess)*contrasts.con(conind).vector];
                                    end
                                elseif ischar(contrasts.con(conind).vector) % contrast string
                                    convec = [convec zeros(size(contrasts.con(conind).vector,1), numcolsinthissess)];
                                    for cr = 1:size(contrasts.con(conind).vector,1)
                                        cs = textscan(contrasts.con(conind).vector(cr,:),'%fx%[^mp]%c%d','Delimiter','|');
                                        if isempty(cs{3}) % simple format
                                            cs = textscan(contrasts.con(conind).vector(cr,:),'%fx%s','Delimiter','|');
                                            cs{3}(1:numel(cs{1}),1) = 'm';
                                            cs{4}(1:numel(cs{1}),1) = int32(1);
                                        end
                                        for e = 1:numel(cs{1})
                                            switch cs{3}(e)
                                                case 'm' % main trail
                                                    EVpttrn = sprintf('Sn(%d) %s*',find(selected_sessions==sess),deblank(cs{2}{e}));
                                                case 'p' % parametric
                                                    EVpttrn = sprintf('Sn(%d) %sx',find(selected_sessions==sess),deblank(cs{2}{e}));
                                            end
                                            ind = find(contains(SPM.xX.name,EVpttrn));
                                            if isempty(ind), ind = find(contains(SPM.xX.name,EVpttrn(1:end-1))); end % covariates
                                            convec(cr,ind(cs{4}(e))) = sessforcon(sess)*cs{1}(e);
                                        end
                                    end
                                else
                                    aas_log(aap,true,'ERROR: Contrast vector must be either string or row vector of numbers.');
                                end
                            else
                                convec = [convec zeros(size(contrasts.con(conind).vector,1), numcolsinthissess)];
                            end
                        end
                        
                        % If subjects have different # of sessions, then
                        % they will be weighted differently in 2nd level
                        % model. So, normalize the contrast by the number
                        % of sesisons that contribute to it [default]
                        
                        if isempty(aas_getsetting(aap,'scalebynumberofsessions')) || aas_getsetting(aap,'scalebynumberofsessions')
                            convec = convec ./ nnz(sessforcon); 
                        end
                        
                    case 'uniquebysession'
                                               
                        totnumcolsbarconstants = size(SPM.xX.X,2) - nsess;
        
                       if ischar(contrasts.con(conind).vector) 
                           
                           % apply string contrast across all sessions as a whole
                           % (events in contrast can be missing in some sessions) 
                           
                            convec = zeros(1,totnumcolsbarconstants);
                                
                            for cr = 1:size(contrasts.con(conind).vector,1)
                                cs = textscan(contrasts.con(conind).vector(cr,:),'%fx%[^mpn]%c%d','Delimiter','|'); 
                                if isempty(cs{3}) % simple format
                                    cs = textscan(contrasts.con(conind).vector(cr,:),'%fx%s','Delimiter','|');
                                    cs{3}(1:numel(cs{1}),1) = 'm';
                                    cs{4}(1:numel(cs{1}),1) = int32(1);
                                end
                                for e = 1:numel(cs{1})
                                    switch cs{3}(e)
                                        case 'm' % main trail
                                            EVpttrn = sprintf(' %s*',deblank(cs{2}{e}));
                                        case 'p' % parametric
                                            EVpttrn = sprintf(' %sx',deblank(cs{2}{e}));
                                        case 'n' % NEW -- contrast can specify nusiance columns using "n" suffix
                                            EVpttrn = sprintf(' %s',deblank(cs{2}{e}));
                                    end
                                    ind = cell_index(SPM.xX.name,EVpttrn);
                                    convec(ind) = cs{1}(e);
                                end
                            end

                       else         
                        
                            if (size(contrasts.con(conind).vector,2) > totnumcolsbarconstants)
                                aas_log(aap,true,sprintf('ERROR: Number of columns in contrast matrix is more than number of columns in model - wanted %d columns, got %d',totnumcolsbarconstants, size(contrasts.con(conind).vector,2)));
                            elseif (size(contrasts.con(conind).vector,2) < totnumcolsbarconstants)
                                convec = [];
                                convec(SPM.xX.iC) = 0;
                                convec(1:numel(contrasts.con(conind).vector)) = contrasts.con(conind).vector;
                                if settings.automatic_movesandmeans
                                    % [AVG] *better* way of specifying the correct columns...
                                    convec_out = zeros(1,totnumcolsbarconstants);
                                    convec_out(SPM.xX.iC) = convec;                                    
                                    convec = convec_out;
                                end
                            else
                                convec=contrasts.con(conind).vector;
                            end
                        
                        
                       end                      
                        
                        
                    otherwise
                        aas_log(aap,true,sprintf('Unknown format %s specified for contrast %d',contrasts.con(conind).format,ccount));
                        
                end
                
                cons{ccount} = [convec zeros(size(convec,1),nsess)];  % Add final constant terms
                
                % Check not empty
                if (~any(cons{ccount}(:)))
                    aas_log(aap,true,sprintf('Contrast %d has no non-zero values, not permitted.',ccount));
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
                
                % contrast name
                aas_log(aap,false,contrasts.con(conind).name)
                
                % list of its event weights

                % if there are a large number of contrasts (e.g., item analysis), echoing
                % each contrast vector to the command window can take quite awhile...
                % -- let the user opt-out

                if ~isempty(aas_getsetting(aap,'diagnostics.display_convecs')) && aas_getsetting(aap,'diagnostics.display_convecs')
                    for cindex  =  1:max(size(convec_names))
                        aas_log(aap,false,sprintf('\t%s: %d', convec_names{cindex}, convec(SPM.xX.iC(cindex))))
                    end
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
            if isempty(SPM.xCon)
                SPM.xCon = spm_FcUtil('Set', connames{conind}, contype{conind},'c', cons{conind}', SPM.xX.xKXs);
            else
                SPM.xCon(end+1) = spm_FcUtil('Set', connames{conind}, contype{conind},'c', cons{conind}', SPM.xX.xKXs);
            end
        end
        SPM = spm_contrasts(SPM);      
        
        % save diagnostic summary images
        diag(aap, subj, SPM);
  
        % save a plaintext list of the contrasts
        writetable(cell2table(convec_names'), fullfile(aas_getsubjpath(aap,subj),'CONTRAST_NAMES.txt'),'WriteVariableNames',0)
  
        % Describe outputs
        
        %  updated spm
        aap = aas_desc_outputs(aap,subj,'firstlevel_spm',fullfile(anadir,'SPM.mat'));
        
        %  firstlevel_betas (includes related statistical files)
        filters={'con','spmT','spmF'};
        
        for filterind=1:length(filters)
            allbetas=dir(fullfile(anadir,[filters{filterind} '_*.nii']));
            betafns=[];
            for betaind=1:length(allbetas)
                betafns=strvcat(betafns,fullfile(anadir,allbetas(betaind).name));
            end
            allbetas=dir(fullfile(anadir,[filters{filterind} '_*.img']));
            for betaind=1:length(allbetas)
                betafns=strvcat(betafns,fullfile(anadir,allbetas(betaind).name));
            end
            allbetas=dir(fullfile(anadir,[filters{filterind} '_*.hdr']));
            for betaind=1:length(allbetas)
                betafns=strvcat(betafns,fullfile(anadir,allbetas(betaind).name));
            end
            aap=aas_desc_outputs(aap,subj,['firstlevel_' lower(filters{filterind}) 's'],betafns);
        end
        cd (cwd);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
        
end
end



function h = diag(aap,subj,SPM)

% Based on Rik Henson's script

% Note this calculation of efficiency takes the 'filtered and whitened'
% design matrix (X) as it is in SPM.

% SPM gets passed in now, because this will load the input stream SPM file,
% which doesn't (always) have the contrasts.

if nargin < 3 % SPM is not passed (e.g. reporting)
    load(aas_getfiles_bystream(aap,subj,'firstlevel_spm'));
end

% distribution
cons = SPM.xCon; cons = cons([cons.STAT]=='T');
if isempty(aas_getsetting(aap,'diagnostics.histogram')) || aas_getsetting(aap,'diagnostics.histogram')
    for c = cons
        h = img2hist(fullfile(SPM.swd, c.Vspm.fname), [], strrep(c.name,' ',''), 0.1);
        print(h,'-djpeg','-r150', fullfile(aas_getsubjpath(aap,subj), ...
            ['diagnostic_' mfilename '_dist_' strrep_multi(c.name,{' ' ':' '>'},{'' '_' '-'}) '.jpg'])); % "unconventional" characters
        close(h);
    end
end

% efficiency
X = SPM.xX.xKXs.X;
iXX=inv(X'*X);

[~, nameCols] = strtok(SPM.xX.name(SPM.xX.iC),' ');
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

% reduce font size as the number of contrasts to display increases
% -- this isn't perfect, but it works better than using a fixed font size

FS = 16;
if (numel(cons)>5);FS=14;end
if (numel(cons)>20);FS=12;end
if (numel(cons)>50);FS=8;end
if (numel(cons)>100);FS=5;end

winheight = 200*numel(cons);
if (winheight > 720); winheight=720; end

h = figure; set(h, 'Position', [1 1 1280 winheight]); % HD
ht = text(1,1,nameCons,'FontSize',FS,'FontWeight','Bold','interpreter','none');
set(ht,'Unit','normalized');
tSize = get(ht,'Extent'); tWidth = tSize(3);
vt = text(1,1,nameCols,'FontSize',FS,'FontWeight','Bold','interpreter','none');
set(vt,'Unit','normalized');
tSize = get(vt,'Extent'); tHeight = 0.1 + tSize(3);

tHeight = 2*tHeight;

close(h);

h = figure; set(h, 'Position', [1 1 1280 winheight]); % HD
subplot('Position', [tWidth 0.15 0.6-tWidth 0.9-tHeight]); % assume not more then 20 contrast
imagesc(columnsCon);
colormap(vertcat(create_grad([0 0 1],[1 1 1],128),create_grad([1 1 1],[1 0 0],128)));
caxis([-max(abs(columnsCon(:))) max(abs(columnsCon(:)))]);
set(gca, 'YTick', 1:numel(cons), 'YTickLabel',strrep(nameCons,'_','-'),  ...
    'Xtick', 1:length(nameCols), 'XTickLabel',strrep(nameCols,'_','-'))
set(gca, 'XAxisLocation','top');

if (numel(cons)>5); set(gca,'XTickLabelRotation',90);end

set(gca,'FontSize',FS,'FontWeight','Bold');
subplot('Position', [0.6 0.15 0.35 0.9-tHeight]);
set(gca, 'YTick', 1:numel(cons), 'YTickLabel','');
set(gca,'FontSize',FS,'FontWeight','Bold');
hold on;
cmap = colorcube(numel(cons));

if aas_getsetting(aap,'estimateefficiency')
    for conind = 1:numel(cons)
        barh(numel(cons) - conind + 1, log(effic(conind)), 'FaceColor',cmap(conind,:));
    end
    ylim([0.5 numel(cons)+0.5])
    title('Log Efficiency','FontSize',18)
    Xs = xlim;
    efficiencyVals = create_grad(Xs(1),Xs(2),5);
    set(gca, 'Xtick', efficiencyVals, 'XtickLabel', compose("%1.1f", exp(efficiencyVals)), 'FontSize', 18)
end

fname = fullfile(aas_getsubjpath(aap,subj),['diagnostic_' mfilename '.jpg']);
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
