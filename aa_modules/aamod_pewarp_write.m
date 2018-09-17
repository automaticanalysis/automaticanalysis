% AA module - coregistration of EPI to structural with unwarping
% [aap,resp]=aamod_pewarp_write(aap,task,subj,sess)
% Ideally preceded by aamod_coreg_extended, to increase rigid alignment!
%% README
%   1) After running aamod_coregister_unwarp the functional and structural
% images are coregistered, so you can directly feed them into
% normalisation.
%   2) It is assumed that the structural has been brain extracted (with BET,
% BSE, etc.), and ideally bias corrected. The mean EPI is bias corrected
% within the script 
%   3) This has only been tested with SPM8. 
%   4) mex_pewarpcost_regularised is a compiled function that is built for 
% Matlab 7.8. It will probably not work in older versions. If you want to 
% recompile you can use the buildC99ompblas script, but it's probably quite
% a hassle as you'll need a recent version of gcc.
% 
% If you run into any problems, let us know!
% Tool developer: Eelke (eelke.visser@donders.ru.nl)
% Ported to AA: Alejandro (a.vicente.grab@gmail.com)
% 
% If you use this tool, please cite Eelke's work:
% EPI DISTORTION CORRECTION BY CONSTRAINED NONLINEAR COREGISTRATION IMPROVES GROUP FMRI
% E. Visser1,2, S. Qin1,3, and M. P. Zwiers1,2
% 1Donders Institute for Brain, Cognition and Behaviour, Radboud University Nijmegen, Nijmegen, Netherlands, 2Department of Psychiatry, Radboud
% University Nijmegen Medical Centre, Nijmegen, Netherlands, 3Department of Neurology, Radboud University Nijmegen Medical Centre, Nijmegen,
% Netherlands
% Proc. Intl. Soc. Mag. Reson. Med. 18 (2010)
%%
function [aap,resp]=aamod_pewarp_write(aap,task,subj,sess)

resp='';

switch task
    case 'doit'
        
        % Load the parameters
        PEparams = []; order = [];
        load(aas_getfiles_bystream(aap,subj,'PEwarp_params'))
        
        %% Now apply this transformation to all the EPI images
        
        % Locate all the EPIs we want to PEwarp
        EPIimg = aas_getfiles_bystream(aap,subj,sess,'epi');
        
        % For each image, apply the warp of the mean EPI image
        aas_log(aap,false,sprintf('PEwarping images for session: %s', aas_getsessdesc(aap,subj,sess)))
        
        PEimg=[];
        for f = 1:size(EPIimg, 1)
            [pth, fn, ext] = fileparts(EPIimg(f,:));
            write_warped_no_jacobian(PEparams, 1e6, order, ...
                EPIimg(f,:), ...
                fullfile(pth, ['p' fn ext])); % Files out come with a 'p' prefix
            PEimg=strvcat(PEimg,fullfile(pth,['p' fn ext]));
        end
        
        %% Describe the outputs
        aas_desc_outputs(aap,subj,sess,'epi',PEimg);
        
    case 'checkrequirements'
        aas_log(aap,0,'Need to trim or skull strip structural\n' );
end