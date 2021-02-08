% AA module - first level statistics for PPI
% [aap,resp]=aamod_ppi_model(aap,task,subj,sess)
% Limitations:
%   - no parametric modulation
% First-level model based on ppi_spm_batch.m 17 by Guillaume Flandin & Darren Gitelman
% Tibor Auer MRC CBSU Jul 2014

function [aap,resp]=aamod_ppi_prepare(aap,task,subj)

resp='';

switch task
    case 'report'
        
        
    case 'doit'
        
        %% Init
        dat = load(aas_getfiles_bystream(aap,'subject',subj,'firstlevel_spm')); SPM = dat.SPM;

        for sess = aap.acq_details.selected_sessions
            fnames = {};
            SPM.pwd = aas_getsesspath(aap,subj,sess); % save folder
            VOIs = cellstr(aas_getfiles_bystream(aap,'session',[subj sess],'vois'));
             
            %% Create PPI
            for p = aas_getsetting(aap,'PPI')                
                fnVOI = VOIs{contains(spm_file(VOIs,'basename'),p.voiname) & strcmp(spm_file(VOIs,'ext'),'mat')};
                
                conspec = p.contrastspec;
                if isempty(regexp(conspec,'[+-][0-9\.]*x', 'once')) % contrast name
                    contrasts = aas_getsetting(aas_setcurrenttask(aap,aas_getsourcestage(aap,'aamod_firstlevel_contrasts','firstlevel_spm')),'contrasts');
                    % - locate contrast specification - no check is needed because it has been already processed
                    contrasts = contrasts(strcmp({contrasts.subject},aas_getsubjname(aap,subj))).con;
                    conspec = contrasts(strcmp({contrasts.name},conspec)).vector;
                end
                
                ppicon = zeros(0,3);
                for ev = strsplit(conspec,'|')
                    e = strsplit(ev{1},'x'); w = str2double(e{1}); e = e{2};
                    ppicon = vertcat(ppicon,[find(strcmp(horzcat(SPM.Sess(sess).U.name),e)) 1 w]);
                end
                
                spm_peb_ppi(SPM,'ppi',fnVOI,ppicon,p.name,1);
                fnames = vertcat(fnames,cellstr(fullfile(aas_getsesspath(aap,subj,sess),['PPI_' p.name '.mat'])));
                spm_print(spm_file(fnames{end},'ext','jpg','prefix',['diagnostic_' mfilename]),spm_figure('GetWin','PPI'),'jpg')
                spm_figure('Close','PPI');
            end
            %% Describe outputs
            aap = aas_desc_outputs(aap,'session',[subj sess],'ppi',fnames);
        end
        
    case 'checkrequirements'

end