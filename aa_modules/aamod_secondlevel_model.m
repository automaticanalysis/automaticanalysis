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
            
            so.img.vol = spm_vol(fullfile(aap.directory_conventions.spmdir,aap.directory_conventions.T1template));
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
            spm_print(['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_groupmasksummary.jpg'],so.figure,'jpg');

            % design
            fSPM = aas_getfiles_bystream(aap,aap.tasklist.currenttask.outputstreams.stream{1}); fSPM = deblank(fSPM(1,:)); %all models are the same (number of inputs may vary)
            load(fSPM);            
            spm_DesRep('DesOrth',SPM.xX);
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
        nsub=length(aap.acq_details.subjects);
        aas_log(aap,false,sprintf('%d subjects',nsub));
        % New option to allow suffix to output file in extraparameters
        if (isfield(aap.tasklist.currenttask.extraparameters,'stats_suffix'))
            stats_suffix=aap.tasklist.currenttask.extraparameters.stats_suffix;
        else
            stats_suffix=[];
        end;
        
        % And make analysis directory
        rfxrootdir = fullfile(aap.acq_details.root,[aap.directory_conventions.rfx stats_suffix]);
        if ~exist(rfxrootdir,'file'); mkdir(aap.acq_details.root,[aap.directory_conventions.rfx stats_suffix]);end
        cd(rfxrootdir);
        
        % Now check all subjects have same number of contrasts and same
        %   contrast names at first level
        clear flSPM
        clear flSPMfn;
        for m=1:nsub
            flSPMfn{m}=aas_getfiles_bystream(aap,m,'firstlevel_spm');
            confiles{m}=aas_getfiles_bystream(aap,m,'firstlevel_cons');
            SPMtemp=load(flSPMfn{m});
            flSPM{m}.SPM.xCon=SPMtemp.SPM.xCon;
%             if (m~=1)
%                 if (length(flSPM{m}.SPM.xCon)~=length(flSPM{1}.SPM.xCon))
%                     aas_log(aap,1,sprintf('Number of contrasts in first level analysis for subject %d different from subject 1. They must be the same for aamod_model_secondlevel to work\n',m));
%                     for n=1:length(flSPM(m).SPM.xCon)
%                         if (flSPM{m}.SPM.xCon(n).name~=flSPM{1}.SPM.xCon(n).name);
%                             aas_log(aap,1,sprintf('Names of contrasts at first level different. Contrast %d has name %s for subject %d but %s for subject 1. They must be the same for aamod_model_secondlevel to work\n',n,flSPM{m}.SPM.xCon(n).name,m,flSPM{1}.xCon(n).name));
%                         end;
%                     end;
%                 end;
%             end;
        end;
        %                phs = 1; conname='UF_S'
        allSPMs = {};
        allbetas = {};
        for n=1:length(flSPM{1}.SPM.xCon)
            
            if strcmp(flSPM{1}.SPM.xCon(n).STAT, 'T')
                
                conname=flSPM{1}.SPM.xCon(n).name;
                % take out characters that don't go well in filenames...
                conname = char(regexp(conname,'[a-zA-Z0-9_-]','match'))';
                rfxdir = fullfile(rfxrootdir,conname);
                if exist(rfxdir)~=7; mkdir(rfxrootdir,conname);end
                cd(rfxdir);
                
                clear SPM
                
                %-Assemble SPM structure
                %=======================================================================
                
                SPM.nscan = nsub;
                SPM.swd = rfxdir;
                
                for s=1:nsub
                    foundit=false;
                    for fileind=1:size(confiles{s},1);
                        [pth nme ext]=fileparts(confiles{s}(fileind,:));
                        m = regexp(nme,sprintf('con_%04d',n));
                        if m
                            foundit=true;
                            break;
                        end;
                    end;
                    if (~foundit)
                        aas_log(aap,true,sprintf('Contrast %d not found in subject %s',n,aap.acq_details.subjects(s).subjname));
                    end;
                    SPM.xY.P{s}   = confiles{s}(fileind,:);
                    SPM.xY.VY(s)   = spm_vol(SPM.xY.P{s});
                end
                
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
                
                % Estimate parameters
                %===========================================================================
                % avoid overwrite dialog
                prevmask = spm_select('List',SPM.swd,'^mask\..{3}$');
                if ~isempty(prevmask)
                    for f = 1 : size(prevmask, 1)
                        spm_unlink(fullfile(SPM.swd, prevmask(f,:)));
                    end
                end
                
                SPM = spm_spm(SPM);
                
                
                % Output streams
                %  secondlevel_spm
                allSPMs{end+1} = fullfile(rfxdir,'SPM.mat');
                
                %  secondlevel_betas (includes related statistical files)
                allbetas{end+1} = char(...
                    spm_select('FPList',rfxdir,'^beta_.*'),...
                    spm_select('FPList',rfxdir,'^ResMS.*'),...
                    spm_select('FPList',rfxdir,'^RPV.*'),...
                    spm_select('FPList',rfxdir,'^mask.*')...
                    );
            end
        end
        
        %% Describe outputs
        aap=aas_desc_outputs(aap,'secondlevel_spm',char(allSPMs));
        aap=aas_desc_outputs(aap,'secondlevel_betas',char(allbetas));
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
end