% AA module
% First-level model Adam Hampshire MRC CBU Cambridge Feb 2006
% Modified for aa by Rhodri Cusack Mar 2006-2011
% Additions by Rik Henson Mar 2011
% Tibor Auer MRC CBU Cambridge 2012-2013

function [aap,resp]=aamod_firstlevel_contrasts(aap,task,subj)

resp='';

switch task
    case 'domain'
        resp='subject';   % this module needs to be run once per subject
        
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
        
        cwd=pwd;
        % get the subdirectories in the main directory
        subj_dir = aas_getsubjpath(aap,subj);
        
        % Maintained for backwards compatibility- better now now put
        % module-specific value in
        % aap.directory_conventions.stats_singlesubj
        if (isfield(aap.tasklist.currenttask.extraparameters,'stats_suffix'))
            stats_suffix=aap.tasklist.currenttask.extraparameters.stats_suffix;
        else
            stats_suffix=[];
        end;
        anadir = fullfile(subj_dir,[aap.directory_conventions.stats_singlesubj stats_suffix]);
        
        % Now set up contrasts...
        SPM=load(aas_getfiles_bystream(aap,subj,'firstlevel_spm'));
        SPM=SPM.SPM;
        SPM.swd=anadir;
        
        % Load up contrasts from task settings
        [fle subjname ext]=fileparts(subj_dir);
        contrasts_set=find(strcmp({aap.tasklist.currenttask.settings.contrasts.subject},subjname));
        if (isempty(contrasts_set))
            % Try for wildcard
            contrasts_set=find(strcmp({aap.tasklist.currenttask.settings.contrasts.subject},'*'));
            if (isempty(contrasts_set))
                aas_log(aap,true,'Can''t find declaration of what contrasts to use - insert this in a local copy of aamod_firstlevel_contrasts.xml or put into user script');
            end;
        end
        
        contrasts=aap.tasklist.currenttask.settings.contrasts(contrasts_set);
        
        for conind=1:length(contrasts.con)
            
            switch(contrasts.con(conind).format)
                case {'singlesession','sameforallsessions'}
                    if (strcmp(contrasts.con(conind).format,'singlesession'))
                        sessforcon=[strcmp({aap.acq_details.sessions.name},contrasts.con(conind).session)];
                    else
                        % [AVG] To make the selected sessions work...
                        sessforcon=zeros(1,length(aap.acq_details.sessions));
                        for sess=aap.acq_details.selected_sessions
                            sessforcon(sess) = 1;
                        end
                        %sessforcon=ones(1,length(SPM.Sess));
                    end;
                    convec=[];
                    sessnuminspm=1;
                    for sess=aap.acq_details.selected_sessions
                        numcolsinthissess=length(SPM.Sess(sessnuminspm).col);
                        if (sessforcon(sess))
                            if (size(contrasts.con(conind).vector,2) > numcolsinthissess)
                                aas_log(aap,true,sprintf('Number of columns in contrast matrix for session %d is more than number of columns in model for this session - wanted %d columns, got ',sess,numcolsinthissess)); disp(contrasts.con(conind).vector);
                            elseif (size(contrasts.con(conind).vector,2) < numcolsinthissess)
                                convec = [convec contrasts.con(conind).vector zeros(size(contrasts.con(conind).vector,1),numcolsinthissess-size(contrasts.con(conind).vector,2))];
                                aas_log(aap,false,sprintf('Warning: Number of columns in contrast matrix for session %d is less than number of columns in model for this session - wanted %d columns, so padding to ',sess,numcolsinthissess)); disp(convec);
                            else
                                convec=[convec contrasts.con(conind).vector];
                            end
                        else
                            convec=[convec zeros(size(contrasts.con(conind).vector,1),numcolsinthissess)];
                        end;
                        sessnuminspm=sessnuminspm+1;
                    end;
                case 'uniquebysession'
                    totnumcolsbarconstants = size(SPM.xX.X,2) - length(aap.acq_details.selected_sessions);
                    if (size(contrasts.con(conind).vector,2) > totnumcolsbarconstants)
                        aas_log(aap,true,sprintf('Number of columns in contrast matrix for session %d is more than number of columns in model (bar constants) - wanted %d columns, got ',totnumcolsbarconstants)); disp(contrasts.con(conind).vector);
                    elseif (size(contrasts.con(conind).vector,2) < totnumcolsbarconstants)
                        convec = [contrasts.con(conind).vector zeros(size(contrasts.con(conind).vector,1),totnumcolsbarconstants-size(contrasts.con(conind).vector,2))];
                        if (contrasts.automatic_movesandmeans)
                            convec_out=[];
                            ind=1;
                            sessnuminspm=1;
                            for sess=aap.acq_details.selected_sessions
                                numcolsinthissess_withoutmoves=length(SPM.Sess(sessnuminspm).col)-6;
                                newind=ind+numcolsinthissess_withoutmoves;
                                convec_out=[convec_out convec(:,ind:(newind-1)) zeros(size(convec,1),6)];
                                ind=newind;
                                sessnuminspm=sessnuminspm+1;
                            end;
                            convec=convec_out;
                        end;
                        if (size(convec,2) < totnumcolsbarconstants)
                            aas_log(aap,false,sprintf('Warning: Number of columns in contrast matrix for ''uniquebysession'' option is less than number columns in model (bar constants) = %d, so padding to ',totnumcolsbarconstants)); disp(convec);
                        end;
                    else
                        convec=contrasts.con(conind).vector;
                    end
                otherwise
                    aas_log(aap,true,sprintf('Unknown format %s specified for contrast %d',contrasts.con(conind).format,conind));
            end;
            
            cons{conind} = [convec zeros(size(convec,1),length(aap.acq_details.selected_sessions))];  % Add final constant terms
            
            % Check not empty
            if (~any(cons{conind}(:)))
                aas_log(aap,true,sprintf('Contrast %d has no non-zero values, not permitted.',contrasts_set(conind)));
            end;
            
            % Get or make automatic contrast name
            if (isfield(contrasts.con(conind),'name') && ~isempty(contrasts.con(conind).name))
                cname{conind}=contrasts.con(conind).name;
            else
                cname{conind}=sprintf('Con%d',conind);
            end
            
            % Allow F tests
            if (isfield(contrasts.con(conind),'type') && isempty(contrasts.con(conind).type))
                type{conind}='T';
            else
                type{conind}=contrasts.con(conind).type;
            end;
        end;
        
        % Make the con images
        SPM.xCon =[];
        for conind = 1:size(cons,2)
            if isempty(SPM.xCon)
                SPM.xCon = spm_FcUtil('Set',cname{conind},type{conind},'c',cons{conind}',SPM.xX.xKXs);
            else
                SPM.xCon(end+1) = spm_FcUtil('Set',cname{conind},type{conind},'c',cons{conind}',SPM.xX.xKXs);
            end
        end
        spm_contrasts(SPM);
        
		% Efficiency based on Rik Henson's script [TA]
        efficiency(aap,subj);
        
        % Describe outputs
        %  updated spm
        aap=aas_desc_outputs(aap,subj,'firstlevel_spm',fullfile(anadir,'SPM.mat'));
        
        %  firstlevel_betas (includes related statistical files)
        filters={'con','spmT','spmF'};
        
        for filterind=1:length(filters)
            allbetas=dir(fullfile(anadir,[filters{filterind} '_*']));
            betafns=[];
            for betaind=1:length(allbetas);
                betafns=strvcat(betafns,fullfile(anadir,allbetas(betaind).name));
            end;
            aap=aas_desc_outputs(aap,subj,['firstlevel_' lower(filters{filterind}) 's'],betafns);
        end;
        cd (cwd);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
end

function f = efficiency(aap,subj)
% Based on Rik Henson's script

% Note this calculation of efficiency takes the 'filtered and whitened'
% design matrix (X) as it is in SPM.
load(aas_getfiles_bystream(aap,subj,aap.tasklist.currenttask.outputstreams.stream{1}));
X = SPM.xX.xKXs.X;
iXX=inv(X'*X);

Sessions = {aap.acq_details.sessions.name};
sCM = aap.tasklist.currenttask.settings.contrasts(subj+1).con;
% CNames = {sCM.name};
f = figure; set(f,'Position',[0 0 1200 400]); 
a = subplot(1,numel(aap.acq_details.selected_sessions)+1,numel(aap.acq_details.selected_sessions)+1);
set(a,'YTickLabel',[]); set(a,'FontSize',16,'FontWeight','Bold'); hold on;
cmap = colorcube(numel(sCM));

CText{1,1} = sprintf('C');
for s = 1:numel(aap.acq_details.selected_sessions)
    sEV = aap.tasksettings.aamod_firstlevel_model.model(1+(subj-1)*numel(aap.acq_details.selected_sessions)+s).event;
    lCon(s) = numel(sEV) + aap.tasksettings.aamod_firstlevel_model.includemovementpars*6;
    EVs = {sEV.name};
    for e = 1:numel(EVs)
        CText = horzcat(CText,EVs{e});
    end
    if aap.tasksettings.aamod_firstlevel_model.includemovementpars
        CText = horzcat(CText,{'x' 'y' 'z' 'r' 'p' 'j'});
    end
end
for s = 1:numel(aap.acq_details.selected_sessions)
    CText = horzcat(CText,sprintf('S%d',s));
end

for c = 1:numel(sCM)
    % Create Contrast
    CM = zeros(1,size(X,2));
    nSess = find(strcmp(sCM(c).session,Sessions)); 
    if numel(sCM(c).vector) <= lCon(nSess) % single session contrast
        CM(sum(lCon(1:nSess-1))+1:sum(lCon(1:nSess))) = sCM(c).vector;
    else % multi session contrast
        CM(1:numel(sCM(c).vector)) = sCM(c).vector;
    end
    
    % Print Contrast
    CText{c+1,1} = sprintf('%02d',c);
    for cc = 1:numel(CM)
        if CM(cc) > 0
            CText{c+1,1+cc} = sprintf(' %1.1f',CM(cc));
        elseif CM(cc) < 0
            CText{c+1,1+cc} = sprintf('%1.1f',CM(cc));
        else
            CText{c+1,1+cc} = sprintf(' %d',CM(cc));
        end        
    end
    
    % Normalize Contrast
    CM = CM / max(sum(CM(CM>0)),sum(CM(CM<0)));
    
    e(c)=trace(CM*iXX*CM')^-1;    
    barh(numel(sCM)-c,e(c),'FaceColor',cmap(c,:));
end

xl = xlim;
yl = ylim;
ylim([-1 yl(2)])
x0 = -numel(aap.acq_details.selected_sessions)*(xl(2)+1)-0.5; 
dx = abs(x0/(size(CText,2)+1));
for y = 1:size(CText,1)
    for x = 1:size(CText,2)
        text(x0+x*dx,numel(sCM)-y+1,CText{y,x},'FontSize',8);
    end
    if y == 1
        text(x0+(x+1)*dx,numel(sCM)-y+1,'Efficiency','FontSize',8);
    end
end
print(f,'-djpeg','-r150',fullfile(aas_getsubjpath(aap,subj),'diagnostic_aamod_firstlevel_contrast'));
close(f);
end


