function [ aap,resp ] = aamod_rescale_epi(aap, task, subject_index, session_index)
%
% aamod_rescale_epi - standardize an epi volume (voxelwise, across time)
% 
% Change History
%
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
  
 		if (~isempty(aap.tasklist.currenttask.settings.scale))
            
            scale = aap.tasklist.currenttask.settings.scale;
            
        else
            
            % mode1000 scaling
           
            header = spm_vol(input_fname);
            V = spm_read_vols(header);
            temp = V(V>0);
            scale = 1000/mode(temp(:));
        
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