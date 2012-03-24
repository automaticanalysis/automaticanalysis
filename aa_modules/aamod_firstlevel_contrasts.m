% AA module
% First-level model Adam Hampshire MRC CBU Cambridge Feb 2006
% Modified for aa by Rhodri Cusack Mar 2006-2011
% Additions by Rik Henson Mar 2011

function [aap,resp]=aamod_firstlevel_contrasts(aap,task,i)

resp='';

switch task
    case 'domain'
        resp='subject';   % this module needs to be run once per subject
        
    case 'description'
        resp='SPM5 contrasts';
        
    case 'summary'
        subjpath=aas_getsubjpath(i);
        resp=sprintf('Contrasts %s\n',subjpath);
        
    case 'report'
        
    case 'doit'
        
        cwd=pwd;
        % get the subdirectories in the main directory
        subj_dir = aas_getsubjpath(aap,i);
        
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
        SPM=load(aas_getfiles_bystream(aap,i,'firstlevel_spm'));
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
            if length(SPM.xCon)==0
                SPM.xCon = spm_FcUtil('Set',cname{conind},type{conind},'c',cons{conind}',SPM.xX.xKXs);
            else
                SPM.xCon(end+1) = spm_FcUtil('Set',cname{conind},type{conind},'c',cons{conind}',SPM.xX.xKXs);
            end
        end
        spm_contrasts(SPM);
        
        % Describe outputs
        %  updated spm
        aap=aas_desc_outputs(aap,i,'firstlevel_spm',fullfile(anadir,'SPM.mat'));
        
        %  firstlevel_betas (includes related statistical files)
        filters={'con','spmT','spmF'};
        
        for filterind=1:length(filters)
            allbetas=dir(fullfile(anadir,[filters{filterind} '_*']));
            betafns=[];
            for betaind=1:length(allbetas);
                betafns=strvcat(betafns,fullfile(anadir,allbetas(betaind).name));
            end;
            aap=aas_desc_outputs(aap,i,['firstlevel_' lower(filters{filterind}) 's'],betafns);
        end;
        cd (cwd);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;



