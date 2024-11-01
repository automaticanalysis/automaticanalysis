% AA module - first level statistics
% **********************************************************************
% You should no longer need to change this module - you may just
% modify the .xml or model in your user script
% **********************************************************************
% Based on original by FIL London and Adam Hampshire MRC CBU Cambridge Feb 2006
% Modified for aa by Rhodri Cusack MRC CBU Mar 2006-Aug 2007
% Thanks to Rik Henson for various suggestions (modified) [AVG & TA]

function [aap,resp]=aamod_firstlevel_model(aap,task,subj)

resp='';

switch task
    
    case 'report' % [TA]
        if ~exist(fullfile(aas_getsubjpath(aap,subj),'diagnostic_aamod_firstlevel_model_design.jpg'),'file')
            load(aas_getfiles_bystream(aap,subj,aap.tasklist.currenttask.outputstreams.stream{1}));
            spm_DesRep('DesOrth',SPM.xX);
            saveas(spm_figure('GetWin','Graphics'),fullfile(aas_getsubjpath(aap,subj),'diagnostic_aamod_firstlevel_model_design.jpg'));
            close all;
        end
        fdiag = dir(fullfile(aas_getsubjpath(aap,subj),'diagnostic_*.jpg'));
        for d = 1:numel(fdiag)
            aap = aas_report_add(aap,subj,'<table><tr><td>');
            aap=aas_report_addimage(aap,subj,fullfile(aas_getsubjpath(aap,subj),fdiag(d).name));
            aap = aas_report_add(aap,subj,'</td></tr></table>');
        end

    case 'doit'
        
        % rWLS toolbox init
        
        if (strcmp(aap.tasklist.currenttask.settings.autocorrelation,'wls'))
            [ ~, WLS ] = aas_cache_get(aap,'wls');
            WLS.load;
        end
        
        % Get subject directory
        cwd=pwd;
        
        % We can now have missing sessions per subject, so we're going to use only
        % the sessions that are common to this subject and selected_sessions
        
        [~, sessInds] = aas_getN_bydomain(aap, 'session', subj);
        subjSessionI = intersect(sessInds, aap.acq_details.selected_sessions);
        numSess = numel(subjSessionI);
        
        modname = aap.tasklist.currenttask.name;
        modnameind = regexp(modname, '_\d{5,5}$');
        modindex = str2num(modname(modnameind+1:end));
        
        % Add PPI if exist
        
        for sess = subjSessionI
            
            if aas_stream_has_contents(aap,'session',[subj,sess],'ppi')
                
                % if the user has defined ppidef_* covariates
                % explicitly, add them to model
                
                [defCov, ind] = aas_getsetting(aap,'modelC','session',[subj sess]);
                if isstruct(defCov) && any(cellfun(@(x) ~isempty(regexp(x,'^ppidef_.*', 'once')), {defCov.covariate.name}))

                    covToRemove = [];
                    PPIs = cellstr(aas_getfiles_bystream(aap,'session',[subj,sess],'ppi'));
                    for indDefPPI = find(cellfun(@(x) ~isempty(regexp(x,'^ppidef_.*', 'once')), {defCov.covariate.name}))
                        defPPI = defCov.covariate(indDefPPI);
                        covToRemove(end+1) = indDefPPI;
                        fnPPI = PPIs{contains(spm_file(PPIs,'basename'),defPPI.vector')};
                        dat = load(fnPPI); PPI = dat.PPI;                        
                        if ~any(strcmp({aap.tasksettings.(modname(1:modnameind-1))(modindex).modelC(ind).covariate.name},['Phys_' PPI.xY.name]))
                            aap = aas_addcovariate(aap,modname,...
                                aas_getsubjname(aap,subj),aas_getsessname(aap,sess),...
                                ['Phys_' PPI.xY.name],PPI.Y,0,1);
                        end
                        aap = aas_addcovariate(aap,modname,...
                            aas_getsubjname(aap,subj),aas_getsessname(aap,sess),...
                            strrep(defPPI.name,'ppidef_',''),PPI.ppi,0,1);
                    end                 
                    aap.tasksettings.(modname(1:modnameind-1))(modindex).modelC(ind).covariate(covToRemove) = []; % remove original definition
              
                else

                    % if the user has not defined ppidef_* covariates explicity,
                    % just add the 3 regressors from the ppi stream. This
                    % makes basic PPI modeling easier to use while still
                    % providing flexible PPI modeling options to those
                    % who need them.
               
                    temp = load(aas_getfiles_bystream(aap,subj,sess,'ppi'));

                    aap = aas_addcovariate(aap,modname,...
                        basename(aas_getsubjpath(aap,subj)),aap.acq_details.sessions(sess).name,...
                        'PPI',temp.PPI.ppi,0,1);
                    aap = aas_addcovariate(aap,modname,...
                        basename(aas_getsubjpath(aap,subj)),aap.acq_details.sessions(sess).name,...
                        'PSYCH',temp.PPI.P,0,1);
                    aap = aas_addcovariate(aap,modname,...
                        basename(aas_getsubjpath(aap,subj)),aap.acq_details.sessions(sess).name,...
                        'PHYS',temp.PPI.Y,0,1);

                end % if ~isempty(defCov{1})...       
            end % if aas_stream_has_contents(aap,'session',[subj,sess],'ppi')...
        end % for sess = subjSessionI...
    
        % update current sesstings
        aap.tasklist.currenttask.settings.modelC = aap.tasksettings.(modname(1:modnameind-1))(modindex).modelC;
        
        % Prepare basic SPM model...
        [SPM, anadir, files, allfiles, model, modelC] = aas_firstlevel_model_prepare(aap, subj);

        % Get all the nuisance regressors...
        [movementRegs, compartmentRegs, physiologicalRegs, spikeRegs, GLMDNregs] = ...
            aas_firstlevel_model_nuisance(aap, subj, files);

        %% Set up CORE model
        
        cols_nuisance=[];
        cols_interest=[];
        currcol=1;

        sessnuminspm=0;

        for sess = 1:numSess
            sessnuminspm=sessnuminspm+1;

            % Settings
            SPM.nscan(sessnuminspm) = size(files{sess},1);
            SPM.xX.K(sessnuminspm).HParam = aap.tasklist.currenttask.settings.highpassfilter;

            % Set up model
            [SPM, cols_interest, cols_nuisance, currcol] = ...
                aas_firstlevel_model_define(aap, sess, sessnuminspm, SPM, model, modelC, ...
                                                             cols_interest, cols_nuisance, currcol, ...
                                                             movementRegs, compartmentRegs, physiologicalRegs, spikeRegs, GLMDNregs);
        end

        cd (anadir)

        %%%%%%%%%%%%%%%%%%%
        %% DESIGN MATRIX %%
        %%%%%%%%%%%%%%%%%%%
        
        SPM.xY.P = allfiles;
        if aap.tasklist.currenttask.settings.firstlevelmasking && aap.tasklist.currenttask.settings.firstlevelmasking < 1
            spm_get_defaults('mask.thresh',aap.tasklist.currenttask.settings.firstlevelmasking);
        end
        
        if (strcmp(aap.tasklist.currenttask.settings.autocorrelation,'wls'))
            SPMdes = spm_rwls_fmri_spm_ui(SPM);
        else
            SPMdes = spm_fmri_spm_ui(SPM);
        end

        SPMdes.xX.X = double(SPMdes.xX.X);
        
        % DIAGNOSTIC
        subjname = aas_prepare_diagnostic(aap, subj);
        try
            saveas(1, fullfile(aap.acq_details.root, 'diagnostics', ...
                                                [mfilename '__' subjname '.fig']));
        catch
        end

        % now check real covariates and nuisance variables are
        % specified correctly
        SPMdes.xX.iG = cols_nuisance;
        SPMdes.xX.iC = cols_interest;

        % Turn off masking if requested
        if ~aap.tasklist.currenttask.settings.firstlevelmasking
            SPMdes.xM=-inf(size(SPMdes.xX.X,1),1);
        end
        
        % Correct epmty model for sphericity check
        if isempty([SPMdes.Sess.U])
            SPMdes.xX.W  = sparse(eye(size(SPMdes.xY.P,1)));
            SPMdes.xVi.V = sparse(eye(size(SPMdes.xY.P,1)));            
        end
        
        %%%%%%%%%%%%%%%%%%%
        %% ESTIMATE MODEL%%
        %%%%%%%%%%%%%%%%%%%
        
        % avoid overwrite dialog
        prevmask = spm_select('List',SPMdes.swd,'^mask\..{3}$');
        if ~isempty(prevmask)
            for ind=1:size(prevmask,1)
                spm_unlink(fullfile(SPMdes.swd, prevmask(ind,:)));
            end
        end

        % add explicit mask if one is specified
        % (SPM will intersect this with any implicit masking...)
        
        if aas_stream_has_contents(aap, subj,'explicitmask')
           explicitmaskfname = aas_getfiles_bystream(aap, subj,'explicitmask');
           aas_log(aap,false,sprintf('Adding explicit mask %s', explicitmaskfname));
           SPMdes.xM.VM = spm_vol(explicitmaskfname);
        end

        if (strcmp(aap.tasklist.currenttask.settings.autocorrelation,'wls'))
            SPMest = spm_rwls_spm(SPMdes);
        else
            SPMest = spm_spm(SPMdes);
        end
        
        % Saving Residuals
        if isfield(aap.tasklist.currenttask.settings,'writeresiduals') && ~isempty(aap.tasklist.currenttask.settings.writeresiduals)
            aas_log(aap,false,'Writing residuals...');
            VRes = spm_write_residuals(SPMest,aap.tasklist.currenttask.settings.writeresiduals);
            for s = 1:numel(SPM.nscan)
                sesspath = aas_getsesspath(aap,subj,subjSessionI(s));
                fres = char({VRes(sum(SPM.nscan(1:s-1))+1:sum(SPM.nscan(1:s))).fname}');
                for f = 1:size(fres,1)
                    movefile(fres(f,:),sesspath);
                end
                residuals{subjSessionI(s)} = horzcat(repmat([sesspath filesep],[SPM.nscan(s) 1]),fres);
                if aap.options.NIFTI4D
                    spm_file_merge(residuals{subjSessionI(s)},fullfile(sesspath,sprintf('Res-%04d.nii',subjSessionI(s))),0,SPM.xY.RT);
                    fres = cellstr(residuals{subjSessionI(s)});
                    residuals{subjSessionI(s)} = fullfile(sesspath,sprintf('Res-%04d.nii',subjSessionI(s)));                
                    delete(fres{:});
                end                
            end
            aas_log(aap,false,'\bDone.');
        end        
        
        %% Describe outputs
        
        cd (cwd);

        %  firstlevel_spm
        aap=aas_desc_outputs(aap,subj,'firstlevel_spm',fullfile(anadir,'SPM.mat'));

        %  firstlevel_betas (includes related statistical files)
        allbetas=dir(fullfile(anadir,'beta_*'));
        allbetas=vertcat(allbetas,...
            dir(fullfile(anadir,'ResMS.*')),...
            dir(fullfile(anadir,'RPV.*')));
        if aap.tasklist.currenttask.settings.firstlevelmasking
            allbetas=vertcat(allbetas,...
                dir(fullfile(anadir,'mask.*')));
            mask = dir(fullfile(anadir,'mask.*'));
            mask=strcat(repmat([anadir filesep],[numel(mask) 1]),char({mask.name}));
            aap=aas_desc_outputs(aap,subj,'firstlevel_brainmask',mask);            
            % save diagnostic images for QA montage
            save_three_ortho_jpgs(mask,fullfile(aas_getsubjpath(aap,subj),'diagnostic_BRAINMASK'));
        end
        betafns=strcat(repmat([anadir filesep],[numel(allbetas) 1]),char({allbetas.name}));
        aap=aas_desc_outputs(aap,subj,'firstlevel_betas',betafns);
        
        if isfield(aap.tasklist.currenttask.settings,'writeresiduals') && ~isempty(aap.tasklist.currenttask.settings.writeresiduals)
            inputstreams = aas_getstreams(aap,'input');        
            residualstreamname = inputstreams{1};
            for s = subjSessionI
                aap=aas_desc_outputs(aap,subj,s,residualstreamname,residuals{s});
            end
        end

        %% DIAGNOSTICS
        
        %  document rwls results using spm_rwls_resstats 

        if (strcmp(aap.tasklist.currenttask.settings.autocorrelation,'wls'))
            % sometimes rwls looks for the movement params using the wrong          
            % filename, so pass in an expliclit fname if we can
            if (aas_stream_has_contents(aap, subj, 'realignment_parameter'))
                for sindex = 1:numSess
                    moveparfname{sindex} = aas_getfiles_bystream(aap, subj, sindex,'realignment_parameter'); 
                end
                spm_rwls_resstats(SPMest,[],moveparfname);
            else
                spm_rwls_resstats(SPMest);
            end          
            
            h = spm_figure('GetWin', 'Graphics');
            set(h,'Renderer','opengl');
            % the following is a workaround for font rescaling weirdness
            set(findall(h,'Type','text'),'FontSize', 10);
            set(findall(h,'Type','text'),'FontUnits','normalized');
            print(h,'-djpeg','-r150', fullfile(aas_getsubjpath(aap,subj), 'diagnostic_RWLS.jpg'));
            spm_figure('Close','Graphics');
            spm_figure('Close','Interactive'); % rWLS toolbox leaves the interactive window up...

        end
 
        % its convenient to save the design matrix now rather than having to run reporting...

        spm_DesRep('DesMtx',SPMdes.xX,[],SPMdes.xsDes);
        h = spm_figure('GetWin', 'Graphics');
        set(h,'Renderer','opengl');
        % the following is a workaround for font rescaling weirdness
        set(findall(h,'Type','text'),'FontSize', 10);
        set(findall(h,'Type','text'),'FontUnits','normalized');
        print(h,'-djpeg','-r150', fullfile(aas_getsubjpath(aap,subj), 'DESIGN_MATRIX.jpg'));
        spm_figure('Close','Graphics');

        % document columns of design matrix as plaintext (SPM skips labels on the plot)

        writetable(cell2table(SPMest.xX.name'), fullfile(aas_getsubjpath(aap,subj),'REGRESSORS.txt'));

        % document event onsets for verification
         
        fid = fopen(fullfile(aas_getsubjpath(aap,subj),'ONSETS.txt'),'w');

        if (fid > 0) 
            for sessindex = 1:numel(SPMest.Sess)
                fprintf(fid,sprintf('\n\t--- SESSION %0d ---\n',sessindex));
                regressors = SPMest.Sess(sessindex).U;
                for rindex = 1:numel(regressors)
                    fprintf(fid,sprintf('\n\t%s\n\n', regressors(rindex).name{1}));
                    onsets = regressors(rindex).ons;
                    for oindex = 1:length(onsets)
                        fprintf(fid,'\t%g\n', onsets(oindex));
                    end
                end
            end
            fclose(fid);
        end
        
        if ~isempty(SPMdes.xX.iC) % model is not empty
            h = firstlevelmodelStats(anadir, [], spm_select('FPList',anadir,'^mask.*'));
            print(h.regs,'-djpeg','-r150', fullfile(aas_getsubjpath(aap,subj), 'diagnostic_aamod_firstlevel_model_regs.jpg'));
            print(h.betas,'-djpeg','-r150', fullfile(aas_getsubjpath(aap,subj), 'diagnostic_aamod_firstlevel_model_betas.jpg'));
            
            close(h.regs)
            close(h.betas)
        end
        
        % unload toolboxes?
        
        if (strcmp(aap.tasklist.currenttask.settings.autocorrelation,'wls'))
            WLS.unload;
        end     
        
    case 'checkrequirements'

        %% rWLS
        
        if (strcmp(aap.tasklist.currenttask.settings.autocorrelation,'wls'))
            if ~aas_cache_get(aap,'wls'), aas_log(aap,true,'rWLS toolbox not found'); end
        end              
        
        %% Add PPI preparations (if needed)
        [~, sessInds] = aas_getN_bydomain(aap, 'session', subj);
        subjSessionI = intersect(sessInds, aap.acq_details.selected_sessions);
        modname = aap.tasklist.currenttask.name;
        modnameind = regexp(modname, '_\d{5,5}$');
        modindex = str2num(modname(modnameind+1:end));
        
        for sess = subjSessionI
            if aas_stream_has_contents(aap,'ppi')
                [defCov, ind] = aas_getsetting(aap,'modelC','session',[subj sess]);
                
                    if isstruct(defCov) && any(cellfun(@(x) ~isempty(regexp(x,'^ppidef_.*', 'once')), {defCov.covariate.name}))
                                      
                    ppiprepstage = aap.tasklist.main.module(aas_getsourcestage(aap,'aamod_ppi_prepare','ppi'));                        
                    ppis = aap.tasksettings.aamod_ppi_prepare(ppiprepstage.index).PPI;

                    for indDefPPI = find(cellfun(@(x) ~isempty(regexp(x,'^ppidef_.*', 'once')), {defCov.covariate.name}))
                        defPPI = defCov.covariate(indDefPPI).vector;
                        
                        if ischar(defPPI) % already processed
                            continue
                        end
                        
                        if any(strcmp({ppis.name},defPPI.name)) % existing definition
                            if ~strcmp(ppis(strcmp({ppis.name},defPPI.name)).voiname,defPPI.voiname) || ~strcmp(ppis(strcmp({ppis.name},defPPI.name)).contrastspec,defPPI.contrastspec)
                                aas_log(aap,true,['ERROR: PPI ' defPPI.name ' with a different definition already exists']);
                            end
                        else
                            ppis(end+1) = defPPI;
                        end
                        % update covariate
                        cov = defCov.covariate(indDefPPI);
                        cov.vector = defPPI.name';
                        aap.tasksettings.(modname(1:modnameind-1))(modindex).modelC(ind).covariate(indDefPPI) = cov;
                        aap.aap_beforeuserchanges.tasksettings.(modname(1:modnameind-1))(modindex).modelC(ind).covariate(indDefPPI) = cov;
                        aap.internal.aap_initial.tasksettings.(modname(1:modnameind-1))(modindex).modelC(ind).covariate(indDefPPI) = cov;
                        aap.internal.aap_initial.aap_beforeuserchanges.tasksettings.(modname(1:modnameind-1))(modindex).modelC(ind).covariate(indDefPPI) = cov;
                    end
                    aap.tasksettings.aamod_ppi_prepare(ppiprepstage.index).PPI = ppis;
                    aap.aap_beforeuserchanges.tasksettings.aamod_ppi_prepare(ppiprepstage.index).PPI = ppis;
                    aap.internal.aap_initial.tasksettings.aamod_ppi_prepare(ppiprepstage.index).PPI = ppis;
                    aap.internal.aap_initial.aap_beforeuserchanges.tasksettings.aamod_ppi_prepare(ppiprepstage.index).PPI = ppis;                    
                end
            end
        end
        
        %% Adjust outstream?
                
        % if we don't write residuals, we need to delete the residual output stream
        % the residual streamname is the same as the input streamname (e.g, "epi" in => "epi" out as residuals)
        % however, we can't assume the name = "epi" because the stream is renameable
        % instead, we assume the name is the first input stream
        
        inputstreams = aas_getstreams(aap,'input');        
        residualstreamname = inputstreams{1};
        if isempty(aas_getsetting(aap,'writeresiduals')) && any(strcmp(aas_getstreams(aap,'output'),residualstreamname))
            aap = aas_renamestream(aap,aap.tasklist.currenttask.name,residualstreamname,[],'output');
            aas_log(aap,false,sprintf('REMOVED: %s output stream: %s', aap.tasklist.currenttask.name,residualstreamname));
        end

        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
        
end
