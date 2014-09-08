% AA module - second level statistics
% Only runs if all contrasts present in same order in all subjects at first
% level. If so, makes model with basic t-test for each of contrasts.
% Second-level model from Rik Henson
% Modified for aa by Rhodri Cusack May 2006
% Tibor Auer MRC CBU Cambridge 2012-2013

function [aap,resp]=aamod_secondlevel_model(aap,task,i)

resp='';

switch task
    case 'domain'
        resp='study';   % this module needs to be run once per study
        
    case 'description'
        resp='SPM5 second level (RFX) model';
        
    case 'summary'
        subjpath=aas_getsubjpath(i);
        resp=sprintf('Second level model %s\n',subjpath);
        
    case 'report'
        if ~exist(fullfile(aas_getstudypath(aap),['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_design.jpg']),'file')
            load(aas_getfiles_bystream(aap,aap.tasklist.currenttask.outputstreams.stream{1}));
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
        
        
        global defaults
        global UFp
        UFp=0.001;
        
        defaults.modality='FMRI'; % Some problems with the Results otherwise?
        
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
        
        for n=1:length(flSPM{1}.SPM.xCon)
            
            if strcmp(flSPM{1}.SPM.xCon(n).STAT, 'T')
                
                conname=flSPM{1}.SPM.xCon(n).name;
                % take out characters that don't go well in filenames...
                conname(conname==':')=[];
                conname(conname==' ')=[];
                conname(conname=='/')=[];
                conname(conname=='\')=[];
                rfxdir = fullfile(rfxrootdir,conname);
                if exist(rfxdir)~=7; mkdir(rfxrootdir,conname);end
                cd(rfxdir);
                
                clear SPM
                
                %-Assemble SPM structure
                %=======================================================================
                
                SPM.nscan = nsub;
                
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
                        aas_log(aap,true,sprintf('Contrast %d not found in subject %s',n,aap.acq_details.subjects(s).mriname));
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
                spm_unlink(fullfile('.', 'mask.img')); % avoid overwrite dialog
                SPM = spm_spm(SPM);
                
                
                % Output streams
                % Describe outputs
                %  secondlevel_spm
                aap=aas_desc_outputs(aap,'secondlevel_spm',fullfile(rfxdir,'SPM.mat'));
                
                %  secondlevel_betas (includes related statistical files)
                allbetas=dir(fullfile(rfxdir,'beta_*'));
                betafns=[];
                for betaind=1:length(allbetas);
                    betafns=strvcat(betafns,fullfile(rfxdir,allbetas(betaind).name));
                end;
                otherfiles={'mask.hdr','mask.img','ResMS.hdr','ResMS.img','RPV.hdr','RPV.img'};
                for otherind=1:length(otherfiles)
                    betafns=strvcat(betafns,fullfile(rfxdir,otherfiles{otherind}));
                end;
                aap=aas_desc_outputs(aap,'secondlevel_betas',betafns);
            end
        end;
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
end