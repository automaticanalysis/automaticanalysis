
function [aap,resp] = aamod_vois_extract(aap, task, subj, sess)
%
% AA module - extract VOIs
%
% based on ppi_spm_batch.m 17 by Guillaume Flandin & Darren Gitelman
% Tibor Auer MRC CBSU Jul 2014
%
% See Chapter 36 (Psychophysiological Interactions) in SPM12 manual
%

resp='';

switch task
    
    case 'domain'
        resp='session';   % this module needs to be run once per subject 
            
    case 'doit'
        
        %% Init
        % Directory
        subj_dir = aas_getsubjpath(aap,subj);
        anadir = fullfile(subj_dir,aap.directory_conventions.stats_singlesubj);
        sessdir = aas_getsesspath(aap,subj,sess);
        
        savedir = pwd;
        cd(anadir);
        delete('VOI*.mat');

        VOIs0 = aap.tasklist.currenttask.settings.VOI;
        fSPM = aas_getfiles_bystream(aap, subj,'firstlevel_spm');
        load(fSPM);
        spm_jobman('initcfg');

        %% Run

        nVOI = 0;
        runF = false;
        
        for v = 1:numel(VOIs0)
            
            VOI = VOIs0(v);
            
            if (VOI.Ic == -1)
               runF = true;
               VOI.Ic = numel(SPM.xCon)+1;
            end
            
            VOI.Sess = sess;
 
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
                        VOI.name = basename(fROI);
                        roi = spm_read_vols(spm_vol(fROI));
                        Vcon = SPM.xCon(VOI.lmax).Vspm;
                        Z = spm_read_vols(Vcon);
                        Z = Z.*roi;
                        lmax_v = find3D(Z==max(Z(Z~=0)))';
                        lmax_mm = Vcon.mat(1:3,:)*[lmax_v; ones(1,size(lmax_v,2))];
                        VOI.xyz = lmax_mm;
                    end
                    nVOI = nVOI + 1;
                    VOIs(nVOI) = VOI;
                end
                
            else
                
                % sanity checks
                
                if isempty(VOI.def)
                    aas_log(aap, true, 'You must define VOI.def as sphere, box, or mask');
                end
                
                % xyz is unused if the VOI is defined as a mask
                % but VOI.xyz can't be empty (and user may have left it
                % empty if passing a mask) 
                
                if (strcmp(VOI.def,'mask'))
                    VOI.xyz = [0 0 0];
                end
                
                if isempty(VOI.spec)
                    aas_log(aap, true, 'You must define VOI.spec as the ROI radius or a mask filename');
                end
                
                if isempty(VOI.lmax)
                    aas_log(aap, true, 'You must supply a valid contrast number in VOI.lmax');
                end
                
                if (isempty(VOI.xyz) || (length(VOI.xyz) ~= 3))
                    aas_log(aap, true, 'You must define VOI.xyz as the ROI center (mm)');
                end

                if isempty(VOI.name)
                    VOI.name = 'unnamed_VOI';
                end
                
                VOI.xyz = VOI.xyz(:);  % spm_regions expects column vector
                nVOI = nVOI + 1;
                VOIs(nVOI) = VOI;
                
            end
                
        end % loop over VOIs0
        
        if runF % generate an F-contrast testing for all effects of interest (and only effects of interest)
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
            
            iSPM = SPM;
            iSPM.title = sprintf('Extracting %s: %s',VOI.name,SPM.xCon(VOI.lmax).name);
            iSPM.Ic = VOI.lmax; % this sets the contrast displayed in the graphics window but is otherwise ignored
            iSPM.n = 1;
            iSPM.Im = [];
            iSPM.pm = [];
            iSPM.Ex = [];
            iSPM.u = 1; % no effective threshold
            iSPM.k = 0;
            iSPM.thresDesc = 'none';
            [hReg,xSPM,SPM] = spm_results_ui('Setup',iSPM);
            
            spm_regions(xSPM,SPM,hReg,VOI);
            
            vfname = sprintf('VOI_%s_%i.mat',VOI.name,VOI.Sess);
            movefile(vfname,sessdir);
            voifn{v} = fullfile(sessdir,vfname);
            
            % save a diagnostic image
            h = spm_figure('GetWin', 'Graphics');
            set(h,'Renderer','opengl');
            % the following is a workaround for font rescaling weirdness
            set(findall(h,'Type','text'),'FontSize', 10);
            set(findall(h,'Type','text'),'FontUnits','normalized');
            print(h,'-djpeg','-r150', fullfile(sessdir, sprintf('VOI_%s_%i.jpg',VOI.name,VOI.Sess)));
            
        end       
        
        % clean up
        
        spm_figure('Close','Interactive');
        spm_figure('Close','Graphics');
        
        %% Describe outputs

        aap=aas_desc_outputs(aap,subj,sess,'vois',char(voifn));
        
        cd(savedir);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;