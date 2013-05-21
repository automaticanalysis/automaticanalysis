% AA module - coregister structural to mean EPI
% Coregistration of structural to mean EPI output by realignment
% Does not require skull stripping any more
% Modified for sparse imaging since prefix for mean is different
% i=subject num
% Rhodri Cusack MRC CBU 2004-6 based on original by Matthew Brett
% 
% Major changes Aug 2010: removed support for central store of structrual
% images. This code was very long in tooth, and unloved.
function [aap,resp] = aamod_coreg_noss(aap, task, subjInd)

resp='';

switch task
    case 'doit'
        global defaults;
        flags = defaults.coreg;
                
        % get mean EPI stream
        PG = aas_getimages_bystream(aap, subjInd,1,'meanepi');
        VG = spm_vol(PG);
        
        % Get path to structural for this subject
        inStream = aap.tasklist.currenttask.inputstreams.stream{1};
        structImg = aas_getfiles_bystream(aap, subjInd, inStream);                
        VF = spm_vol(structImg);

        % do coregistration
        x  = spm_coreg(VG, VF,flags.estimate);
        
        M  = inv(spm_matrix(x));
          
        spm_get_space(structImg, M*spm_get_space(structImg));
       
        aap = aas_desc_outputs(aap, subjInd, inStream, structImg);

        % Save graphical output - this will now be done by report task
        try
            figure(spm_figure('FindWin', 'Graphics'));
        catch
            figure(1);
        end
        print('-djpeg','-r75',fullfile(aas_getsubjpath(aap, subjInd),'diagnostic_aamod_coreg'));
        
    case 'checkrequirements'
        
end
