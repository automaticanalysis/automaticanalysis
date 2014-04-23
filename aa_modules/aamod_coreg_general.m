% AA module - Coregister and reslice images (general)
% i=subject num
% jt (03/July/2012) based loosely on Rhodri's aamod_coreg_noss
% 

function [aap,resp]=aamod_coreg_general(aap,task,i)

resp='';

switch task
    case 'doit'

        % Load default coregistration parameters:
        flags = aap.spm.defaults.coreg;
        flags.write.which = [1 0];    % don't reslice first image or mean
        flags.write.mean = false;     % don't compute/write mean
        % Experimental 
        %  + smaller separation iteration; - tolerances, + interpolation:
        flags.estimate.sep = [4 2 1];
        flags.estimate.tol = .10*flags.estimate.tol;
        flags.write.interp = 7;
            
        % Get target image:
        targetimfn = aas_getimages_bystream(aap,i,[],aap.tasklist.currenttask.settings.inputstreams.stream{1});
        Vtarget = spm_vol(targetimfn);
        
        % Get image to coregister ('source'):
        sourceimfn = aas_getimages_bystream(aap,i,[],aap.tasklist.currenttask.settings.inputstreams.stream{2});
        Vsource = spm_vol(sourceimfn);
        
        % Estimate coregistration parameters:
        x = spm_coreg(Vsource,Vtarget,flags.estimate);
        
        % Apply coregistration parameters to source file:
        spm_get_space(sourceimfn,spm_matrix(x)*Vsource.mat);
                        
        % Reslice:
        spm_reslice(char({targetimfn;sourceimfn}),flags.write);
        
        % Describe outputs:
        [pth fstem fext] = fileparts(sourceimfn);
        reslicedimfn = fullfile(pth,['r' fstem fext]);
        aas_desc_outputs(aap,i,sprintf('%s',aap.tasklist.currenttask.settings.outputstreams.stream),reslicedimfn);
        
        % (Graphical output should be saved by report task)
        
    case 'checkrequirements'

end
