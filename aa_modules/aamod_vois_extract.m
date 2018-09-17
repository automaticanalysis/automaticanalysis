% AA module - extract VOIs
% [aap,resp]=aamod_vois_extract(aap,task,subj,sess)
% based on ppi_spm_batch.m 17 by Guillaume Flandin & Darren Gitelman
% Tibor Auer MRC CBSU Jul 2014

function [aap,resp]=aamod_vois_extract(aap,task,subj,sess)

resp='';

switch task
    case 'domain'
        resp='session';   % this module needs to be run once per subject case 'report'
        
    case 'doit'
        
        %% Init
        % Directory
        subj_dir = aas_getsubjpath(aap,subj);
        anadir = fullfile(subj_dir,aap.directory_conventions.stats_singlesubj);
        sessdir = aas_getsesspath(aap,subj,sess);
        
        % Settings
        cd(anadir);
        delete('VOI*.mat');

        VOIs0 = aap.tasklist.currenttask.settings.VOI;
        fSPM = aas_getfiles_bystream(aap, subj,'firstlevel_spm');
        load(fSPM);
        spm_jobman('initcfg');

        %% Run
        % Setup
        nVOI = 0;
        runF = false;
        for v = 1:numel(VOIs0)
            VOI = VOIs0(v);
            if (VOI.Ic == -1)
               runF = true;
               VOI.Ic = numel(SPM.xCon)+1;
            end
            VOI.Sess = sess; % single session only
            if isempty(VOI.name), VOI.name = VOI.xyz; end % default

            if isempty(VOI.lmax) && ~strcmp(VOI.def,'mask') % specified
                nVOI = nVOI + 1;
                VOIs(nVOI) = VOI;
            else
                if aas_stream_has_contents(aap,'rois')
                    VOI0 = VOI;
                    ROIs = aas_getfiles_bystream(aap,subj,'rois');
                    for r = 1:size(ROIs,1)
                        VOI = VOI0;
                        fROI = deblank(ROIs(r,:));
                        if strcmp(VOI.def,'mask')
                            VOI.name = basename(fROI);
                            VOI.spec = spm_vol(fROI);
                        else % lmax
                            VOI.defIc = VOI.lmax;
                            VOI.name = basename(fROI);
                            roi = spm_read_vols(spm_vol(fROI));
                            
                            Vcon = SPM.xCon(VOI.defIc).Vspm;
                            Z = spm_read_vols(Vcon);
                            Z = Z.*roi;
                            lmax_v = find3D(Z==max(Z(Z~=0)))';
                            lmax_mm = Vcon.mat(1:3,:)*[lmax_v; ones(1,size(lmax_v,2))];
                            
                            VOI.xyz = lmax_mm;
                        end
                        nVOI = nVOI + 1;
                        VOIs(nVOI) = VOI;
                    end
                end
            end
        end
        
        if runF % generate an F-contrast testing for all effects of interest
            clear jobs
            jobs{1}.stats{1}.con.spmmat = cellstr(fSPM);
            jobs{1}.stats{1}.con.consess{1}.fcon.name = 'Effects of interest';
            fcont = [eye(numel(SPM.Sess(1).U)) zeros(numel(SPM.Sess(1).U),size(SPM.Sess(1).C.C,2)+1)];
            for i=1:size(fcont,1)
                jobs{1}.stats{1}.con.consess{1}.fcon.convec{1,i} = fcont(i,:);
            end
            spm_jobman('run',jobs);
            load(fSPM);
        end
        
        for v = 1:numel(VOIs)
            VOI = VOIs(v);
            
            if ~isfield(VOI,'defIc')
                VOI.defIc = VOI.Ic;
            end
            
            % Extract
            iSPM = SPM;
            iSPM.title = sprintf('Extracting %s: %s',...
                VOI.name,...
                SPM.xCon(VOI.defIc).name);
            iSPM.Ic = VOI.defIc;
            iSPM.n = 1;
            iSPM.Im = [];
            iSPM.pm = [];
            iSPM.Ex = [];
            iSPM.u = 1; % no effective threshold
            iSPM.k = 0;
            iSPM.thresDesc = 'none';
            [hReg,xSPM,SPM] = spm_results_ui('Setup',iSPM);
            
            [Y,xY]  = spm_regions(xSPM,SPM,hReg,VOI);
            
            vfname = sprintf('VOI_%s_%i.mat',VOI.name,VOI.Sess);
            movefile(vfname,sessdir);
            voifn{v} = fullfile(sessdir,vfname);
        end
        
        %% Describe outputs
        %  firstlevel_spm
        aap=aas_desc_outputs(aap,subj,sess,'vois',char(voifn));
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;