% AA module - first level statistics
% **********************************************************************
% You should no longer need to change this module - you may just
% modify the .xml or model in your user script
% **********************************************************************
% Based on original by FIL London and Adam Hampshire MRC CBU Cambridge Feb 2006
% Modified for aa by Rhodri Cusack MRC CBU Mar 2006-Aug 2007
% Thanks to Rik Henson for various suggestions

function [aap,resp]=aamod_GLMdenoise(aap,task,subj)

resp='';

switch task
    case 'report'
        
    case 'doit'
        epifn= aas_getfiles_bystream(aap,subj,1, 'epi');
        
        % Get mask
        Mimg = aas_getfiles_bystream(aap,subj, 'epiBETmask');
        Mimg = Mimg(1,:); % Only first one
        
        % Check orientation of mask and epi the same
        V=spm_vol(char(epifn(1,:),Mimg(1,:)));
        if ~spm_check_orientations(V)
            aas_log(aap,true,'Image orientation mismatch. You may need to reslice your meanepi before epi bet reslice - do you have aamod_norm_write_meanepi in your tasklist?');
        end;
        
        % Load mask
        M = logical(spm_read_vols(spm_vol(Mimg)));
        
        % Session Split
        session_split = aap.tasklist.currenttask.settings.session_split;
        if isempty(session_split)
            session_split{1} = aap.acq_details.selected_sessions;
        end
        % Stimulus Duration in seconds...
        stimdur = aap.tasklist.currenttask.settings.stimdur;
        % Options for GLMdenoise
        opt = aap.tasklist.currenttask.settings.opt;
        optFN = fieldnames(opt);
        for o = 1:length(optFN)
            if isempty(opt.(optFN{o}))
                opt = rmfield(opt, optFN{o});
            end
        end
        opt.brainmask = M;
        
        for z = 1:length(session_split)
            aapSPM = aap;
            aapSPM.acq_details.selected_sessions = session_split{z};
            % Prepare basic SPM model...
            [SPM, anadir, files, allfiles, model, modelC] = aas_firstlevel_model_prepare(aapSPM, subj);
            TR = SPM.xY.RT;
            
            gd_data = cell(1, length(session_split{z}));
            for s = 1:length(session_split{z});
                sess = session_split{z}(s);
                aas_log(aap, 0,  sprintf('Loading gd_data and model of sess %d', sess));
                
                % Get gd_data
                V = spm_vol(files{sess}(1,:));
                %gd_data{s} = single(nan(V.dim(1), V.dim(2), V.dim(3), size(files{sess},1)));
                gd_data{s} = single(nan(sum(M(:)), size(files{sess},1)));
                for f = 1:size(files{sess},1)
                    %gd_data{s}(:,:,:,f) = spm_read_vols(spm_vol(files{sess}(f,:)));
                    Y = spm_read_vols(spm_vol(files{sess}(f,:)));
                    gd_data{s}(:,f) = Y(M);
                end
            end
            
%             memtoc
            
            switch aap.tasklist.currenttask.settings.GDmode
                case ''
                    aas_log(aap, 1, 'You must specify the GDmode parameter, which sets how we use GLMdenoise')
                case 'onsets'
                    hrfmodel = 'assume';
                    hrfknobs = [];
                    
                    clear SPM
                    
                    if isempty(stimdur)
                        aas_log(aap, 1, 'You should specify the stimulus duration (in seconds) in your recipe');
                    end
                    
                    gd_design = cell(1, length(session_split{z}));
                    for s = 1:length(session_split{z});
                        sess = session_split{z}(s);
                        
                        % Set up model
                        gd_ons = cell(1, length(model{sess}.event));
                        for e = 1:length(model{sess}.event)
                            ons = (model{sess}.event(e).ons - 1) * TR; % in seconds & -1 to be in same coordinate system as GLMdenoise
                            dur = ceil(model{sess}.event(e).dur * TR ./ stimdur); % in seconds / stimulus duration
                            dur(dur == 0) = 1;
                            
                            gd_ons{e} = [];
                            for o = 1:length(ons)
                                for d = 1:dur(o);
                                    gd_ons{e} = [gd_ons{e}; ons(o) + (d - 1) * stimdur];
                                end
                            end
                        end
                        
                        gd_design{s} = gd_ons;
                    end
                    
                case 'SPMdesign'
                    hrfmodel = 'assume';
                    hrfknobs = 1;
                    
                    if isempty(stimdur)
                        aas_log(aap, 0, 'Stimdur is set to TR');
                        stimdur = TR;
                    end
                    
                    % Get all the nuisance regressors...
                    [movementRegs, compartmentRegs, physiologicalRegs, spikeRegs] = ...
                        aas_firstlevel_model_nuisance(aapSPM, subj, files);
                    
                    %% Set up CORE model
                    
                    cols_nuisance=[];
                    cols_interest=[];
                    currcol=1;
                    sessnuminspm=0;
                    
                    for sess = session_split{z}
                        sessnuminspm=sessnuminspm+1;
                        
                        % Settings
                        SPM.nscan(sessnuminspm) = size(files{sess},1);
                        SPM.xX.K(sessnuminspm).HParam = aap.tasklist.currenttask.settings.highpassfilter;
                        
                        % Set up model
                        [SPM, cols_interest, cols_nuisance, currcol] = ...
                            aas_firstlevel_model_define(aap, sess, sessnuminspm, SPM, model, modelC, ...
                            cols_interest, cols_nuisance, currcol, ...
                            movementRegs, compartmentRegs, physiologicalRegs, spikeRegs);
                    end
                    
                    cd (anadir)
                    
                    %%%%%%%%%%%%%%%%%%%
                    %% DESIGN MATRIX %%
                    %%%%%%%%%%%%%%%%%%%
                    
                    SPM.xY.P = allfiles;
                    SPMdes = spm_fmri_spm_ui(SPM);
                    
                    % DIAGNOSTIC
                    mriname = aas_prepare_diagnostic(aap, subj);
                    try
                        saveas(1, fullfile(aap.acq_details.root, 'diagnostics', ...
                            [mfilename '__' mriname '.fig']));
                    catch
                    end
                    
                    %%%%%%%%%%%%%%%%%%%
                    %% GLMdenoise    %%
                    %%%%%%%%%%%%%%%%%%%
                    
                    chunks = [0 cumsum(SPMdes.nscan)];
                    gd_design = cell(size(SPMdes.Sess));
                    
                    for s = 1:length(SPMdes.Sess);
                        rows = (chunks(s) + 1):chunks(s+1);
                        cols = SPMdes.Sess(s).col;
                        cols = cols(ismember(cols, cols_interest));
                        
                        gd_design{s} = SPMdes.xX.X(rows, cols);
                        
                    end
            end
            
            cd(anadir); % So that figures are printed to the right location
            
            % Try without denoising first...
            optT = opt;
            optT.numpcstotry = 0;
            
            gd_dataORI = gd_data;
            
            [gd_resultsALT, gd_dataALT] = ...
                GLMdenoisedata(gd_design, gd_data, stimdur, TR, hrfmodel, hrfknobs, optT, sprintf('figures%d_ALT', z));
            
            % Try with denoising...
            [gd_results, gd_data] = ...
                GLMdenoisedata(gd_design, gd_data, stimdur, TR, hrfmodel, hrfknobs, opt, sprintf('figures%d', z));
            
            [h, figName] = GLMdenoise_diagnostics(gd_results, gd_data, gd_resultsALT, gd_dataORI);
            for f = 1:length(h);
                print(h(f), fullfile(anadir, sprintf('%s_%02d.eps', figName{f}, z)), '-depsc')
                close(h(f))
            end
            
            %% SAVE DENOISED gd_dataALT TO DISC
            files_denoised = cell(size(files));
            for s = 1:length(session_split{z});
                sess = session_split{z}(s);
                
                for f = 1:size(files{sess},1)
                    V = spm_vol(files{sess}(f,:));
                    Y(M) = gd_data{s}(:,f);
                    
                    % Write out denoised file to different filename!
                    [pth, fn, ext] = fileparts(V.fname);
                    V.fname = fullfile(pth, ['d', fn, ext]);
                    spm_write_vol(V,Y);
                    
                    files_denoised{sess} = strvcat(files_denoised{sess}, V.fname);
                end
                
                aap=aas_desc_outputs(aap,subj,sess, 'epi', files_denoised{sess});
            end
            
            % And save gd_results
%             gdfn=fullfile(anadir,'gd_results.mat');
            gdfn=fullfile(aas_getsubjpath(aap, subj), 'gd_results.mat');
            save(gdfn,'gd_results');
%             aap=aas_desc_outputs(aap,subj,sess,'gd_results',gdfn);
            aap=aas_desc_outputs(aap,subj,'gd_results',gdfn);
            
        end
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end