function [aap,resp]=aamod_secondlevel_randomise_threshold(aap, task)
%
% thresholding for aamod_secondlevel_randomise
%
% CHANGE HISTORY
%
% 06/2021 [MSJ] - cleanup
% 05/2020 [MSJ] - filenaming simplifcation and washu rendering options
% 09/2019 [MSJ] - new
%

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

        corrected_tmap_suffix = aap.tasklist.currenttask.settings.corrected_tmap_suffix;
        renderer = aap.tasklist.currenttask.settings.renderer;
       
        % loop over the FSL files, threshold the tmap using the pmap, and save jpegs
               
        for tindex = 1:size(fsl_tmaps,1)
            
            tmap_fname = deblank(fsl_tmaps(tindex,:));            
            tmap_header = spm_vol(tmap_fname);
            
            % lookup pmap for this tmap according to pmap option
      
            [ results_directory,tname,~ ] = fileparts(tmap_fname);
            
            % aamod_secondlevel_randomise uses double underscore for uniqueness
            temp = split(tname,'__');  % NB: double underscore
             
            if strcmp(pmap_option,'uncorrected')
                pname = [ temp{1} '__' temp{2} '__' temp{3} ];
            else
                pname = [ temp{1} '__' temp{2} '__' pmap_option '_corrp_' temp{3} ];      
            end                     
            
            pmap_fname = fullfile(results_directory,[ pname '.nii' ]);
            
            % its possible the file won't exist if the selected pmap option isn't
            % consistent with options passed to aamod_secondlevel_randomise
            
           if (isempty(dir(pmap_fname)))
                aas_log(aap, true, sprintf('Can''t find pmap %s. Check selected thresholding option against the flags passed to Randomise', pmap_fname));
            end
                       
            P = char(pmap_fname, tmap_fname);
            f = sprintf('(i1>%f).*i2',1-thresh); % careful: FSL p maps are 1-p       
            Q = tmap_header;          

            % for one-sample ttest, there will only be *_tstat1.nii which we
            % would like to rename simply to *_[corrected_tmap_suffix].nii
            % without numbering. For more complicated models, there may be 
            % *_tstat1.nii, *_tstat2.nii, etc. and we need to keep the numbering. 
            % As a compromise, we omit the number on the first tstat and
            % keep the numbering on the others (if present)
            
            if (contains(tname,'tstat1'))
                Q.fname =  fullfile(results_directory, [ strrep(tname, 'tstat1', corrected_tmap_suffix) '.nii']);
            else
                Q.fname =  fullfile(results_directory,[ strrep(tname, 'tstat', corrected_tmap_suffix) '.nii']);
            end
            
            flags = {[],[],[],[]};

            % note i1 = pname and i2 = tname
            
            Q = spm_imcalc(P, Q, f, flags);
            
            if (strcmp(renderer,'classic') || strcmp(renderer,'both'))     
                caption = aap.tasklist.currenttask.settings.description;
                label = sprintf('%s \n(%s < %0.4g) .* %s', caption, pname, thresh, tname);
                [p,n,~] = fileparts(Q.fname);
                savefile_fname = fullfile(p,n);
                overlay_nifti(Q.fname, template_fname, render_fname, savefile_fname, label);
            end
           
            if (strcmp(renderer,'washu') || strcmp(renderer,'both')) 
                cfg = [];
                eval([aap.tasklist.currenttask.settings.render_options ';']);
                warning('off','MATLAB:subscripting:noSubscriptsSpecified');
                washu_surfacerender(Q.fname,[],cfg,Q.fname);        
            end            
 
            threshmaps{tindex} = Q.fname;
 
        end
        
        if (aap.tasklist.currenttask.settings.do_file_cleanup)

            for tindex = 1:size(fsl_tmaps,1)
                tmap_fname = deblank(fsl_tmaps(tindex,:));
                [results_directory,oldname,~] = fileparts(tmap_fname);
                if (contains(oldname,'tstat1'))
                    newname = strrep(oldname,'tstat1','uncorrected');
                else
                    newname = strrep(oldname,'tstat','uncorrected');
                end
                aas_shell(sprintf('mv %s/%s.nii %s/%s.nii', results_directory, oldname, results_directory, newname));
                aas_shell(sprintf('rm -f %s/*_corrp_*.nii', results_directory));
                oldname = strrep(oldname,'tstat1','glm_cope');
                if (exist(fullfile(results_directory,[oldname '.nii']),'file'))
                    newname = strrep(oldname,'_glm_cope','_cope');
                    aas_shell(sprintf('mv %s/%s.nii %s/%s.nii', results_directory, oldname, results_directory, newname));
                    aas_shell(sprintf('rm -f %s/*glm*.nii', results_directory));
                end
             end

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




