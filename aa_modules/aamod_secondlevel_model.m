% AA module - second level statistics
% Only runs if all contrasts present in same order in all subjects at first
% level. If so, makes model with basic t-test for each of contrasts.
% Second-level model from Rik Henson
% Modified for aa by Rhodri Cusack May 2006
% Tibor Auer MRC CBU Cambridge 2012-2013

function [aap,resp]=aamod_secondlevel_model(aap,task)

resp='';

switch task
    case 'report'
        [~, SPMtool] = aas_cache_get(aap,'spm');
        
        if numel(cellstr(spm_select('List',aas_getstudypath(aap),'^diagnostic.*'))) < 2 % two images expected
            % group mask
            for subj = 1:numel(aap.acq_details.subjects)
                mask{subj} = aas_getfiles_bystream(aap,'subject',subj,'firstlevel_brainmask');
            end
            Vm = cell2mat(spm_vol(mask));
            Ym = spm_read_vols(Vm);
            Ys = sum(Ym,4);
            Vs = Vm(1);
            Vs.fname = fullfile(aas_getstudypath(aap),'groupmaps_summary.nii');
            spm_write_vol(Vs,Ys)
            
            so = slover;
            so.figure = spm_figure('GetWin', 'SliceOverlay');
            
            so.img.vol = spm_vol(fullfile(SPMtool.toolPath,aap.directory_conventions.T1template));
            so.img.prop = 1;
            so.transform = 'axial';
            so = fill_defaults(so);
            so.slices = -72:6:108;

            so.img(2).vol = Vs;
            so.img(2).type  = 'truecolour';
            so.img(2).prop  = 0.33;
            so.img(2).cmap  = hot;
            so.img(2).range = [0 max(Ys(:))];
            so.cbar = [so.cbar 2];

            so = paint(so);
            spm_print(fullfile(aas_getstudypath(aap),['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_groupmasksummary.jpg']),so.figure,'jpg');

            % design
            fSPM = aas_getfiles_bystream(aap,aap.tasklist.currenttask.outputstreams.stream{1}); fSPM = deblank(fSPM(1,:)); %all models are the same (number of inputs may vary)
            dat = load(fSPM);
            spm_DesRep('DesOrth',dat.SPM.xX);
            saveas(spm_figure('GetWin','Graphics'),fullfile(aas_getstudypath(aap),['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_design.jpg']));
            close all;
        end
        fdiag = dir(fullfile(aas_getstudypath(aap),'diagnostic_*.jpg'));
        for d = 1:numel(fdiag)
            aap = aas_report_add(aap,[],'<table><tr><td>');
            aap=aas_report_addimage(aap,[],fullfile(aas_getstudypath(aap),fdiag(d).name));
            aap = aas_report_add(aap,[],'</td></tr></table>');
        end
    case 'doit'
        %% Init
        nsub=aas_getN_bydomain(aap,'subject');
        aas_log(aap,false,sprintf('%d subjects',nsub));
        % New option to allow suffix to output file in extraparameters
        if (isfield(aap.tasklist.currenttask.extraparameters,'stats_suffix'))
            stats_suffix=aap.tasklist.currenttask.extraparameters.stats_suffix;
        else
            stats_suffix=[];
        end
        
        % And make analysis directory
        rfxrootdir = fullfile(aap.acq_details.root,[aap.directory_conventions.rfx stats_suffix]);
        if ~exist(rfxrootdir,'file'); mkdir(aap.acq_details.root,[aap.directory_conventions.rfx stats_suffix]);end
        cd(rfxrootdir);
        
        %% Setup model
        SPM.nscan = nsub;
        SPM.xX = struct(	'X',	ones(nsub,1),...
            'iH',1,'iC',zeros(1,0),'iB',zeros(1,0),'iG',zeros(1,0),...
            'name',{{'mean'}},'I',[[1:nsub]' ones(nsub,3)],...
            'sF',{{'obs'  ''  ''  ''}});
        
        SPM.xC = [];
        
        SPM.xGX = struct(...
            'iGXcalc',1,	'sGXcalc','omit',				'rg',[],...
            'iGMsca',9,	'sGMsca','<no grand Mean scaling>',...
            'GM',0,		'gSF',ones(nsub,1),...
            'iGC',	12,	'sGC',	'(redundant: not doing AnCova)',	'gc',[],...
            'iGloNorm',9,	'sGloNorm','<no global normalisation>');
        
        SPM.xVi	= struct('iid',1,'V',speye(nsub));
        
        Mdes 	= struct(	'Analysis_threshold',	{'None (-Inf)'},...
            'Implicit_masking',	{'Yes: NaNs treated as missing'},...
            'Explicit_masking',	{'Yes: SPM2 Brain Mask'});
        
        %SPM.xM	= struct('T',-Inf,'TH',ones(nsub*2,1)*-Inf,...
        %		 'I',1,'VM',spm_vol('/home/rh01/SPM/spm5/apriori/brainmask.nii'),'xs',Mdes);
        
        SPM.xM	= struct('T',-Inf,'TH',ones(nsub*2,1)*-Inf,...
            'I',1,'VM',[],'xs',Mdes);
        
        Pdes 	= {{'1 condition, +0 covariate, +0 block, +0 nuisance'; '1 total, having 1 degrees of freedom'; 'leaving 8 degrees of freedom from 9 images'}};
        
        SPM.xsDes = struct(	'Design',		{'One sample t-test'},...
            'Global_calculation',	{'omit'},...
            'Grand_mean_scaling',	{'<no grand Mean scaling>'},...
            'Global_normalisation',	{'<no global normalisation>'},...
            'Parameters',		Pdes);
        
        %% Add data
        if aas_stream_has_contents(aap,'subject',1,'firstlevel_spm')
            % Now check all subjects have same number of contrasts and same
            %   contrast names at first level
            clear flSPM
            clear flSPMfn;
            for m=1:nsub
                flSPMfn{m}=aas_getfiles_bystream(aap,m,'firstlevel_spm');
                confiles{m}=aas_getfiles_bystream(aap,m,'firstlevel_cons');
                SPMtemp=load(flSPMfn{m});
                flSPM{m}.SPM.xCon=SPMtemp.SPM.xCon;
                if m > 1
                    if numel(flSPM{m}.SPM.xCon) < numel(flSPM{1}.SPM.xCon)
                        aas_log(aap,1,sprintf('Number of contrasts in first level analysis for subject %d smaller than that for subject 1. There MUST be corresponding contrasts for aamod_model_secondlevel to work\n',m));
                        for n=1:numel(flSPM(m).SPM.xCon)
                            if ~strcmp(flSPM{m}.SPM.xCon(n).name,flSPM{1}.SPM.xCon(n).name)
                                aas_log(aap,1,sprintf('Names of contrasts at first level different. Contrast %d has name %s for subject %d but %s for subject 1. They MUST be the same for aamod_model_secondlevel to work\n',n,flSPM{m}.SPM.xCon(n).name,m,flSPM{1}.xCon(n).name));
                            end
                        end
                    end
                end
            end
            %                phs = 1; conname='UF_S'
            flCon = flSPM{1}.SPM.xCon;
            flCon = flCon(strcmp({flCon.STAT},'T'));
            for n=1:length(flCon)
                conname=flCon(n).name;
                % take out characters that don't go well in filenames...
                conname = char(regexp(conname,'[a-zA-Z0-9_-]','match'))';
                rfxdir = fullfile(rfxrootdir,conname);
                aas_makedir(aap, rfxdir);
                cd(rfxdir);
                
                SPM(n) = SPM(1);
                
                %-Assemble SPM structure
                %=======================================================================
                SPM(n).swd = rfxdir;
                
                for s=1:nsub
                    foundit=false;
                    for fileind=1:size(confiles{s},1)
                        if contains(spm_file(confiles{s}(fileind,:),'basename'),sprintf('con_%04d',n))
                            foundit=true;
                            break;
                        end
                    end
                    if (~foundit)
                        aas_log(aap,true,sprintf('Contrast %d not found in subject %s',n,aap.acq_details.subjects(s).subjname));
                    end
                    SPM(n).xY.P{s} = confiles{s}(fileind,:);
                    SPM(n).xY.VY(s) = spm_vol(SPM(n).xY.P{s});
                end
            end
        elseif startsWith(aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name,'aamod_tdt') % CAVE: assumption on module name
            load(aas_getfiles_bystream(aap,'subject',1,'settings'),'cfg');
            meas = cfg.results.output;
            modNames = meas;
            if aas_stream_has_contents(aap,'subject',1,'settings_pairwise')
                for fnpw = cellstr(aas_getfiles_bystream(aap,'subject',1,'settings_pairwise'))'
                    load(fnpw{1},'cfg');
                    modNames = cat(1,modNames,spm_file(meas,'suffix',['_' strjoin(strrep(unique(cfg.files.labelname,'stable'),' ','-'),'_vs_')]));
                end                
            end
            
            for n = 1:numel(modNames)
                rfxdir = fullfile(rfxrootdir,modNames{n});
                aas_makedir(aap, rfxdir);
                SPM(n) = SPM(1);
                SPM(n).swd = rfxdir;
            end
            
            for subj = 1:nsub
                for n = 1:numel(meas)
                    SPM(n).xY.P{subj} = aas_getfiles_bystream(aap,'subject',subj,meas{n});
                    SPM(n).xY.VY(subj) = spm_vol(SPM(n).xY.P{subj});
                end
                if aas_stream_has_contents(aap,'subject',subj,'settings_pairwise')
                    for n = 1:numel(meas)
                        fnpw = cellstr(aas_getfiles_bystream(aap,'subject',subj,[meas{n} '_pairwise']));
                        for pw = 1:numel(fnpw)  
                            SPM(pw*numel(meas)+1).xY.P(subj) = fnpw(pw);
                            SPM(pw*numel(meas)+1).xY.VY(subj) = spm_vol(SPM(pw*numel(meas)+1).xY.P{subj});
                        end
                    end
                end
            end
        else
            aas_log(aap,true,'Unknown use case');
        end
           
        %% Estimate parameters
        allSPMs = cell(1,numel(SPM));
        allbetas = cell(1,numel(SPM));
        for n = 1:numel(SPM)
            % avoid overwrite dialog
            prevmask = spm_select('List',SPM(n).swd,'^mask\..{3}$');
            if ~isempty(prevmask)
                for f = 1 : size(prevmask, 1)
                    spm_unlink(fullfile(SPM(n).swd, prevmask(f,:)));
                end
            end
            
            spm_spm(SPM(n));
            
            % Output streams
            %  secondlevel_spm
            allSPMs{n} = fullfile(SPM(n).swd,'SPM.mat');
            
            %  secondlevel_betas (includes related statistical files)
            allbetas{n} = char(...
                spm_select('FPList',SPM(n).swd,'^beta_.*'),...
                spm_select('FPList',SPM(n).swd,'^ResMS.*'),...
                spm_select('FPList',SPM(n).swd,'^RPV.*'),...
                spm_select('FPList',SPM(n).swd,'^mask.*')...
                );
        end
        
        %% Describe outputs
        aap = aas_desc_outputs(aap,'study',[],'secondlevel_spm',char(allSPMs));
        aap = aas_desc_outputs(aap,'study',[],'secondlevel_betas',char(allbetas));
    case 'checkrequirements'
        if ~aas_cache_get(aap,'spm'), aas_log(aap,true,'SPM is not found'); end
        
        % tdt
        if startsWith(aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name,'aamod_tdt') % CAVE: assumption on module name
            srcmodulename = 'aamod_tdt_decode'; % CAVE: assumption on module name
            src = aas_getstreams(aas_setcurrenttask(aap,aas_getsourcestage(aap,srcmodulename)),'output'); src = setdiff(src,{'settings' 'mask'});
            inp = aas_getstreams(aap,'input');
            if any(strcmp(src,'settings_pairwise'))
                if ~any(strcmp(inp,'settings_pairwise')), aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'append','settings_pairwise','input'); end
                src = setdiff(src,{'settings_pairwise'});
            end
            inp(contains(inp,{'settings' 'settings_pairwise' 'firstlevel_brainmask'})) = [];
            for s = 1:numel(src) 
                if s == 1
                    if strcmp(inp{s},src{s}), continue; end
                    aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'input',src{s},'input');
                else
                    aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'append',src{s},'input');
                end
                aas_log(aap,false,['INFO: ' aap.tasklist.currenttask.name ' input streams: ''' src{s} '''']);
            end
        end
        
end
end