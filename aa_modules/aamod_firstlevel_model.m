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
        % Get subject directory
        cwd=pwd;
        
 		% We can now have missing sessions per subject, so we're going to use only
        % the sessions that are common to this subject and selected_sessions
        [numSess, sessInds] = aas_getN_bydomain(aap, 'session', subj);
        subjSessionI = intersect(sessInds, aap.acq_details.selected_sessions);
        numSess = numel(subjSessionI);        
        modname = aap.tasklist.currenttask.name;
        modnameind = regexp(modname, '_\d{5,5}$');
        modindex = str2num(modname(modnameind+1:end));
        
        % Add PPI if exist
        for sess = subjSessionI
            if aas_stream_has_contents(aap,'session',[subj,sess],'ppi')
                isDefPPI = true;
                [defPPI, ind] = aas_getsetting(aap,'modelC','session',[subj sess]);
                if isempty(defPPI)
                    aas_log(aap,true,sprintf('WARNING: no PPI specification for %s',aas_getsubjdesc(aap,subj)));
                    isDefPPI = false;
                end
                indDefPPI = cellfun(@(x) ~isempty(regexp(x,'^PPI.*', 'once')), {defPPI.covariate.name});
                if ~any(indDefPPI)
                    aas_log(aap,true,sprintf('WARNING: no PPI specification for %s',aas_getsubjdesc(aap,subj))); 
                    isDefPPI = false;
                end
                if sum(indDefPPI) > 1
                    aas_log(aap,true,sprintf('WARNING: more than 1 PPI specification for %s -> the first will be used',aas_getsubjdesc(aap,subj)));                     
                    indDefPPI = find(indDefPPI,1,'first');
                end
                if isDefPPI
                    defPPI = defPPI.covariate(indDefPPI);
                    aap.tasksettings.(modname(1:modnameind-1))(modindex).modelC(ind).covariate(indDefPPI) = []; % remove original definition
                    
                    PPIs = cellstr(aas_getfiles_bystream(aap,'session',[subj,sess],'ppi'));
                    fnPPI = PPIs{contains(spm_file(PPIs,'basename'),defPPI.vector')};
                    dat = load(fnPPI); PPI = dat.PPI;
                    
                    aap = aas_addcovariate(aap,modname,...
                        aas_getsubjname(aap,subj),aas_getsessname(aap,sess),...
                        defPPI.name,PPI.ppi,0,1);
%                     aap = aas_addcovariate(aap,modname,...
%                         aas_getsubjname(aap,subj),aas_getsessname(aap,sess),...
%                         ['Psych_' psych],PPI.P,0,1);
                    aap = aas_addcovariate(aap,modname,...
                        aas_getsubjname(aap,subj),aas_getsessname(aap,sess),...
                        ['Phys_' PPI.xY.name],PPI.Y,0,1);
                end
            end
        end
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
        SPMdes = spm_fmri_spm_ui(SPM);

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
            end;
        end
                
        SPMest = spm_spm(SPMdes);
        
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

        % Describe outputs
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
        end
        betafns=strcat(repmat([anadir filesep],[numel(allbetas) 1]),char({allbetas.name}));
        aap=aas_desc_outputs(aap,subj,'firstlevel_betas',betafns);
        
        if isfield(aap.tasklist.currenttask.settings,'writeresiduals') && ~isempty(aap.tasklist.currenttask.settings.writeresiduals)
            for s = subjSessionI
                aap=aas_desc_outputs(aap,subj,s,'epi',residuals{s});
            end
        end

        %% DIAGNOSTICS...
        if ~isempty(SPMdes.xX.iC) % model is not empty
            h = firstlevelmodelStats(anadir, [], spm_select('FPList',anadir,'^mask.*'));
            print(h.regs,'-djpeg','-r150', fullfile(aas_getsubjpath(aap,subj), 'diagnostic_aamod_firstlevel_model_regs.jpg'));
            print(h.betas,'-djpeg','-r150', fullfile(aas_getsubjpath(aap,subj), 'diagnostic_aamod_firstlevel_model_betas.jpg'));
            
            close(h.regs)
            close(h.betas)
        end
    case 'checkrequirements'
        %% Adjust outstream
        if isempty(aas_getsetting(aap,'writeresiduals')) && any(strcmp(aas_getstreams(aap,'output'),'epi'))
            aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'epi',[],'output');
            aas_log(aap,false,sprintf('REMOVED: %s output stream: epi', aap.tasklist.currenttask.name'));
        end

    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
