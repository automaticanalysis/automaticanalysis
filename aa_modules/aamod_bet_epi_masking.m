% AA module
% Runs EPI masking after BET

function [aap,resp]=aamod_bet_epi_masking(aap,task,subj,sess)

resp='';

switch task
    case 'summary'
        subjpath=aas_getsubjpath(subj);
        resp=sprintf('Align %s\n',subjpath);
        
    case 'report'
        
    case 'doit'
        
        % Get epis
        EPIimg = aas_getfiles_bystream(aap,subj,sess,'epi');
        % Get mask
        Mimg = aas_getfiles_bystream(aap,subj, 'epiBETmask');
        % Only first one
        Mimg = Mimg(1,:);
        
        % Mask epis with mask!
        mask_img(Mimg, EPIimg);
        
        %% DESCRIBE OUTPUTS!
        aap=aas_desc_outputs(aap,subj,sess,'epi', EPIimg);
end
