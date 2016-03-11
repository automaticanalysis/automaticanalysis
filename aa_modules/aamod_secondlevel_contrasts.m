% AA module
% It makes basic f- and t-tests.
% Tibor Auer MRC CBU Cambridge 2012-2013

function [aap,resp]=aamod_secondlevel_contrasts(aap,task)

resp='';


switch task
    case 'report'
        fdiag = dir(fullfile(aas_getstudypath(aap),'diagnostic_*.jpg'));
        for d = 1:numel(fdiag)
            aap = aas_report_add(aap,[],'<table><tr><td>');
            aap=aas_report_addimage(aap,[],fullfile(aas_getstudypath(aap),fdiag(d).name));
            aap = aas_report_add(aap,[],'</td></tr></table>');
        end
        
    case 'doit'
        nsub=length(aap.acq_details.subjects);
        aas_log(aap,false,sprintf('%d subjects',nsub));
        % New option to allow suffix to output file in extraparameters
        if (isfield(aap.tasklist.currenttask.extraparameters,'stats_suffix'))
            stats_suffix=aap.tasklist.currenttask.extraparameters.stats_suffix;
        else
            stats_suffix=[];
        end;
        
        rfxrootdir = fullfile(aap.acq_details.root,[aap.directory_conventions.rfx stats_suffix]);
        cd(rfxrootdir);
        
        %% Now set up contrasts
        fnSPMs = aas_getfiles_bystream(aap,'secondlevel_spm');
        for flc = 1:size(fnSPMs,1)
            fnSPM = deblank(fnSPMs(flc,:));
            loaded=load(fnSPM);
            SPM=loaded.SPM;
            SPM.swd=fileparts(fnSPM);
            
            % group mean F
            SPM.xCon = [];
            for c = aap.tasklist.currenttask.settings.contrasts
                if isempty(SPM.xCon)
                    SPM.xCon = spm_FcUtil('Set', c.name, c.type, 'c', c.vector', SPM.xX.xKXs);
                else
                    SPM.xCon(end+1) = spm_FcUtil('Set', c.name, c.type, 'c', c.vector', SPM.xX.xKXs);
                end
            end
            
            SPM = spm_contrasts(SPM);
            
            % Output streams
            %  secondlevel_spm
            allSPMs{flc} = fnSPM;
            
            %  secondlevel_betas (includes related statistical files)
            allCons{flc} = spm_select('FPList',SPM.swd,'^con_.*');
            allTs{flc} = spm_select('FPList',SPM.swd,'^spmT_.*');
            allFs{flc} = spm_select('FPList',SPM.swd,'^spmF_.*');
            
            % Diagnostics (check distribution of T-values in contrasts)
            cons = SPM.xCon; cons = cons([cons.STAT]=='T');
            for c = cons
                h = img2hist(fullfile(SPM.swd, c.Vspm.fname), [], strrep(c.name,' ',''), 1);
                print(h,'-djpeg','-r150', fullfile(aas_getstudypath(aap), ...
                    ['diagnostic_aamod_firstlevel_contrast_dist_' basename(SPM.swd) '_' strrep(c.name,' ','') '.jpg']));
                close(h);
            end
            
        end
        
        %% Describe outputs
        %  updated spm
        aap = aas_desc_outputs(aap,'secondlevel_spm',char(allSPMs));
        if ~isempty(char(allCons)), aap = aas_desc_outputs(aap,'secondlevel_cons',char(allCons)); end
        if ~isempty(char(allCons)), aap = aas_desc_outputs(aap,'secondlevel_spmts',char(allTs)); end
        if ~isempty(char(allCons)), aap = aas_desc_outputs(aap,'secondlevel_spmfs',char(allFs)); end        
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
        
end
end