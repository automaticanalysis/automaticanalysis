% AA module - second level thresholding for FSL randomise
%
% 09/2019 - [MSJ] new


function [aap,resp]=aamod_secondlevel_randomise_threshold(aap, task)

resp='';

switch task
	
    case 'report'
 		
    case 'doit'
		     
        module_directory = aas_getstudypath(aap);
        
        savedir = pwd;
        cd(module_directory);
         
        template_fname = 'toolbox/OldNorm/T1.nii';
        render_fname = 'rend/render_single_subj.mat';
        
        fsl_tmaps = aas_getfiles_bystream(aap, 'secondlevel_fslts');

        thresh = aap.tasklist.currenttask.settings.threshold;
        pmap_option = aap.tasklist.currenttask.settings.pmap;
        
        threshmaps = cell(size(fsl_tmaps,1),1);
        
        % loop over the FSL files, threshold the tmap using the pmap, and save jpegs
               
        for tindex = 1:size(fsl_tmaps,1)
            
            tmap_fname = deblank(fsl_tmaps(tindex,:));            
            tmap_header = spm_vol(tmap_fname);

            % lookup pmap for this tmap according to pmap option
            
            [ p,tname,~ ] = fileparts(tmap_fname);
            temp = split(tname,'_');
            pname = [ temp{1} '_' temp{2} '_' pmap_option '_corrp_' temp{3} ];
            pmap_fname = fullfile(p,[ pname '.nii' ]);
            
            if (isempty(dir(pmap_fname)))
                aas_log(aap, true, sprintf('Cannot find probability map %s. Exiting...', pmap_fname));
            end
                       
            P = char(pmap_fname, tmap_fname);
            f = sprintf('(i1>%f).*i2',1-thresh); % careful: FSL p maps are 1-p       
            Q = tmap_header;          
            Q.fname =  fullfile(module_directory, sprintf('%s_THRESH.nii',tname));
            flags = {[],[],[],[]};

            % note i1 = pname and i2 = tname
            
            Q = spm_imcalc(P, Q, f, flags);
     
            description = aap.tasklist.currenttask.settings.description;
            label = sprintf('%s \n(%s < %0.4g) .* %s', description, pname, thresh, tname);

            savefile_fname = tname;
            overlay_nifti(Q.fname, template_fname, render_fname, savefile_fname, label);
 
            threshmaps{tindex} = Q.fname;
            
        end
        
        % desc outputs
		
        aap = aas_desc_outputs(aap,'secondlevel_thr', threshmaps);
        
         % restore pwd and we're done!

        cd(savedir);
       
    case 'checkrequirements'
        
    otherwise
        
        aas_log(aap, true, sprintf('Unknown task %s', task));
        
end

end

