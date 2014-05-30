function [aap,resp] = aamod_GLMdenoise(aap, task, subjInd)

resp='';

switch task
    case 'report'
        
    case 'doit'
%        memtic
        
        settings = aap.tasklist.currenttask.settings;
        
        epiFiles = aas_getfiles_bystream(aap, subjInd, 1, 'epi'); % get 1st session epis
        
        % We can now have missing sessions per subject, so we're going to use only
        % the sessions that are common to this subject and selected_sessions
        [numSess, sessInds] = aas_getN_bydomain(aap, 'session', subjInd);
        subjSessionI = intersect(sessInds, aap.acq_details.selected_sessions);
        
        % Get mask, if it's an input
        if  aas_stream_has_contents(aap,'epiBETmask')
            Mimg = aas_getfiles_bystream(aap,subjInd, 'epiBETmask');
            Mimg = Mimg(1,:); % Only first one
            
            % Check orientation of mask and epi the same
            Vs1 = spm_vol(epiFiles(1,:));
            Vm = spm_vol(Mimg(1,:));
            
            if ~spm_check_orientations([Vs1 Vm])
                aas_log(aap,true,'Image orientation mismatch. You may need to reslice your meanepi before epi bet reslice - do you have aamod_norm_write_meanepi in your tasklist?');
            end
            
            M = ~logical(spm_read_vols(spm_vol(Mimg)));
            brainI = find(~M);
        else
            M = 0;
            Vs1 = spm_vol(epiFiles(1,:));
            brainI = 1:prod(Vs1.dim);
        end
        
        opt.brainexclude = M;
        
        % Session Split
        session_split = aap.tasklist.currenttask.settings.session_split;
        if isempty(session_split)
            session_split{1} = subjSessionI;
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

        for split = 1:length(session_split)
            aapSPM = aap;
            aapSPM.acq_details.selected_sessions = session_split{split};
            
            % Prepare basic SPM model...
            [SPM, anadir, files, allfiles, model, modelC] = aas_firstlevel_model_prepare(aapSPM, subjInd);
            TR = SPM.xY.RT;
            
            gd_data = cell(1, length(session_split{split}));
            
            % Load the data for each session split
            for s = 1:length(session_split{split});
                sess = session_split{split}(s);
                aas_log(aap, 0,  sprintf('Loading gd_data and model of sess %d', sess));
                
                % Get gd_data
                Vs = spm_vol(files{sess});
                for v = 1 : length(Vs)
                    
                    Ys = single(spm_read_vols(Vs(v)));
                
                    gd_data{s}(:,v) = Ys(brainI);
                end
            end
            
%            memtoc
            
            switch aap.tasklist.currenttask.settings.GDmode
                
                case ''
                    aas_log(aap, 1, 'You must specify the GDmode parameter, which sets how we use GLMdenoise')
                    
                case 'onsets'
                    
                    warning('GDmode ''onsets'' has not yet been tested extensively.  You might be better off trying GDmode = SPMdesign');
                    
                    hrfmodel = 'assume';
                    hrfknobs = [];
                    
                    clear SPM
                    
                    if isempty(stimdur)
                        aas_log(aap, 1, 'You should specify the stimulus duration (in seconds) in your recipe');
                    end
                    
                    gd_design = cell(1, length(session_split{split}));
                    for s = 1:length(session_split{split});
                        sess = session_split{split}(s);
                        
                        % Set up model
                        gd_ons = cell(1, length(model{sess}.event));
                        %                         for e = 1:length(model{sess}.event)
                        %                             ons = (model{sess}.event(e).ons);% - 1) * TR; % in seconds & -1 to be in same coordinate system as GLMdenoise
                        %                             dur = ceil(model{sess}.event(e).dur * TR ./ stimdur); % in seconds / stimulus duration
                        %                             dur(dur == 0) = 1;
                        %
                        %                             gd_ons{e} = [];
                        %                             for o = 1:length(ons)
                        %                                 for d = 1:dur(o);
                        %                                     gd_ons{e} = [gd_ons{e}; ons(o) + (d - 1) * stimdur];
                        %                                 end
                        %                             end
                        %                         end
                        
                        gd_design{s} = {model{sess}.event.ons};
                    end
                    
                case 'SPMdesign'
                    hrfmodel = 'assume';
                    hrfknobs = 1;
                    
                    if isempty(stimdur)
                        aas_log(aap, 0, 'Stimdur is set to TR');
                        stimdur = TR;
                    end
                    
                    % Get all the nuisance regressors...
                    [movementRegs, compartmentRegs, physiologicalRegs, spikeRegs, GLMDNregs] = ...
                        aas_firstlevel_model_nuisance(aapSPM, subjInd, files);
                    
                    %% Set up CORE model
                    
                    cols_nuisance=[];
                    cols_interest=[];
                    currcol=1;
                    sessnuminspm=0;
                    
                    for sess = session_split{split}
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
                    SPMdes = spm_fmri_spm_ui(SPM);
                    
                    % DIAGNOSTIC
                    mriname = aas_prepare_diagnostic(aap, subjInd);
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
            
            if settings.generateALTresults
                [gd_resultsALT, gd_dataALT] = ...
                    GLMdenoisedata(gd_design, gd_data, stimdur, TR, hrfmodel, hrfknobs, optT, sprintf('figures%d_ALT', split));
            end
            
            % Try with denoising...
            [gd_results, gd_data] = ...
                GLMdenoisedata(gd_design, gd_data, stimdur, TR, hrfmodel, hrfknobs, opt, sprintf('figures%d', split));
            
            gd_results.models = [];
            
            gdFile = fullfile(aas_getsubjpath(aap, subjInd), 'gd_results.mat');
            save(gdFile, 'gd_results', '-v7.3');
            aap = aas_desc_outputs(aap, subjInd, 'gd_results', gdFile);
            
            %             [h, figName] = GLMdenoise_diagnostics(gd_results, gd_data, gd_resultsALT, gd_dataORI);
            %             for f = 1:length(h);
            %                 print(h(f), fullfile(anadir, sprintf('%s_%02d.eps', figName{f}, z)), '-depsc')
            %                 close(h(f))
            %             end
            
            %% SAVE DENOISED gd_dataALT TO DISC
            
            if ismember({'epi'}, settings.outputstreams.stream)
                
                files_denoised = cell(size(files));
                for s = 1:length(session_split{split});
                    sess = session_split{split}(s);
                    
                    for f = 1:size(files{sess},1)
                        V = spm_vol(files{sess}(f,:));

                        % Write out denoised file to different filename!
                        [pth, fn, ext] = fileparts(V.fname);
                        V.fname = fullfile(pth, ['d', fn, ext]);
                        spm_write_vol(V,gd_data{s}(:,:,:,f));
                        
                        files_denoised{sess} = strvcat(files_denoised{sess}, V.fname);
                    end
                    
                    aap=aas_desc_outputs(aap,subjInd,sess, 'epi', files_denoised{sess});
                end
                
                aap=aas_desc_outputs(aap,subj,sess, 'epi', files_denoised{sess});
            end
        end
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end