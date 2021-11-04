function [ aap,resp ] = aamod_rescale_epi(aap, task, subject_index, session_index)
%
% aamod_rescale_epi - standardize an epi volume (voxelwise, across time)
% 
% Change History
%
% summer 2021 [MSJ] - add explicit_mask
% summer 2020 [MSJ] - change to fslmaths
% spring 2020 [MSJ] - new
%

resp='';

switch task
	
    case 'report'
				
    case 'doit'
        					
		input_fname = aas_getimages_bystream(aap, subject_index, session_index, 'epi');		
 		[ ~,name,~ ] = fileparts(input_fname);
        
        output_fname = fullfile(aas_getsesspath(aap, subject_index, session_index),['z' name '.nii']);
  
        scale = aap.tasklist.currenttask.settings.scale;
        
        % adjust raw scale for median and mode scaling
        
        switch aap.tasklist.currenttask.settings.scalingmode
            
            case 'median'

                header = spm_vol(input_fname);
                V = spm_read_vols(header);

                if ~isempty(aap.tasklist.currenttask.settings.explicit_mask)

                    if strcmp(aap.tasklist.currenttask.settings.explicit_mask,'default')
                        mask_fname = spm_select('FPListRec',aap.directory_conventions.spmdir,'brainmask.nii');
                    else
                        mask_fname = aap.tasklist.currenttask.settings.explicit_mask;
                    end
                    mask_header = spm_vol(mask_fname);
                    M = spm_read_vols(mask_header);
                    temp = V(M>0);

                else                                      
                    temp = V(V>0);
                end

                scale = scale/median(round(temp(:)));
           
            case 'mode'
                          
                header = spm_vol(input_fname);
                V = spm_read_vols(header);

                if ~isempty(aap.tasklist.currenttask.settings.explicit_mask)

                    if strcmp(aap.tasklist.currenttask.settings.explicit_mask,'default')
                        mask_fname = spm_select('FPListRec',aap.directory_conventions.spmdir,'brainmask.nii');
                    else
                        mask_fname = aap.tasklist.currenttask.settings.explicit_mask;
                    end
                    mask_header = spm_vol(mask_fname);
                    M = spm_read_vols(mask_header);
                    temp = V(M>0);

                else                                      
                    temp = V(V>0);
                end

                scale = scale/mode(round(temp(:)));
           
            otherwise
        
        end
        
        odt = aap.tasklist.currenttask.settings.odt;
        
        % fslmaths command example for reference:
        %
        %   fslmaths unscaled.nii -mul 0.005 scaled.nii
        %
        
        command = sprintf('fslmaths %s -mul %f %s -odt %s', input_fname, scale, output_fname, odt);
        aas_runfslcommand(aap, command);
 
        aap = aas_desc_outputs(aap, subject_index, session_index, 'epi', output_fname);

		% done!
		        
    case 'checkrequirements'
        
        % should check here if fsl is installed
        
    otherwise
        aas_log(aap, 1, sprintf('%s:Unknown task %s',mfilename, task));
        
end