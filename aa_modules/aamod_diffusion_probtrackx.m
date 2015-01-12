% Get diffusion data streams to apply probtractx.
% Probtrackx repetitively samples from the distributions on voxel-wise 
% principal diffusion directions, each time computing a streamline through
% these local samples to generate a probabilistic streamline or a sample 
% from the distribution on the location of the true streamline. 
% By taking many such samples FDT is able to build up the posterior 
% distribution on the streamline location or the connectivity distribution.

function [aap resp]=aamod_diffusion_probtrackx(aap,task,subjind,sessind,splitind)
global aaworker
resp='';

switch task
    case 'report'
    case 'doit'
        % Launch all of the bedpost jobs
        fslext=aas_getfslext(aap);        
        
        % Rename diffusion_data input
        %diffusion_data=aas_getfiles_bystream(aap,'diffusion_session',[subjind sessind],'diffusion_data');
        %[pth nme ext]=aas_fileparts(diffusion_data);
        %movefile(diffusion_data,fullfile(pth,['data' fslext]));

        % Paths
        dsesspth= aas_getpath_bydomain(aap,'diffusion_session',[subjind sessind]);
        bedpostpath=[dsesspth '.bedpostX'];
        
        psesspth= aas_getpath_bydomain(aap,'diffusion_session_probtrackx',[subjind sessind splitind]);
        
        targets=aas_getfiles_bystream(aap,'diffusion_session',[subjind sessind],'tractography_targets');        
        targetmaskfn=fullfile(psesspth,'target_mask_list.txt');
        fid=fopen(targetmaskfn,'w');
        for targetind=1:size(targets,1)
            fprintf(fid,'%s\n',deblank(targets(targetind,:)));
        end;
    % Remove bet_ from start of bet_nofid_brain_mask etc as bedpost
            % doesn't like this
            betmaskfns=aas_getfiles_bystream(aap,'diffusion_session',[subjind sessind], 'BETmask');
            
            nodif_brain_mask=betmaskfns(1,:);  % top row    
             %[pth nme ext]=aas_fileparts(nodif_brain_mask);  
           %  if strcmp(nme(1:4),'bet_')
              % nme=nme(5:end);
           %  end;
            % nodif_fn=fullfile(pth, [nme ext]);
                          
        fclose(fid);
        
        nstreamlines=round(aap.tasklist.currenttask.settings.totalstreamlines/aap.options.probtrackx.nsplits);
        
        cmd='probtrackx2';
        cmd=[cmd sprintf(' -x %s',aas_getfiles_bystream(aap,'diffusion_session',[subjind sessind],'tractography_seeds'))];
        cmd=[cmd ' -l --onewaycondition --omatrix2'];
        cmd=[cmd sprintf(' -c 0.2 -S 2000 --steplength=0.5 -P %d --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --forcedir --opd --os2t --rseed=%d',nstreamlines,floor(toc))];
        cmd=[cmd sprintf(' --target2=%s',aas_getfiles_bystream(aap,'diffusion_session',[subjind sessind],'dti_FA'))];
        cmd=[cmd sprintf(' -s %s',fullfile(bedpostpath,'merged'))];
        cmd=[cmd sprintf(' -m %s', nodif_brain_mask)] 
        cmd=[cmd sprintf(' --dir=%s',fullfile(psesspth,'logdir'))];
        cmd=[cmd sprintf(' --targetmasks=%s',targetmaskfn)]

        % Do the task
        [s w]=aas_runfslcommand(aap,cmd);
        if s 
            aas_log(aap,true,sprintf('Error executing\n  %s\nof\n%s',cmd,w));
        end     
        
        % Check the output is there
        if ~exist(fullfile(psesspth,'logdir',sprintf('seeds_to_%d_diffusion_space_bin.nii.gz',size(targets,1))),'file')
            aas_log(aap,true,sprintf('Missing output for\n  %s\n',cmd));
        end;
        
        % Write output streams
        fns={};
        for targetind=1:size(targets,1)
            fns{targetind}=fullfile(psesspth,'logdir',sprintf('seeds_to_%d_diffusion_space_bin.nii.gz',targetind));
        end;
        
        
        aap=aas_desc_outputs(aap, 'diffusion_session_probtrackx',[subjind sessind splitind],'seeds_to_diffusion_space',fns);
        
        aap=aas_desc_outputs(aap, 'diffusion_session_probtrackx',[subjind sessind splitind],'fdt_matrix2',{fullfile(psesspth,'logdir','fdt_matrix2.dot'),fullfile(psesspth,'logdir','tract_space_coords_for_fdt_matrix2'),  fullfile(psesspth,'logdir','coords_for_fdt_matrix2')});
        aap=aas_desc_outputs(aap, 'diffusion_session_probtrackx',[subjind sessind splitind],'fdt_paths',fullfile(psesspth,'logdir',['fdt_paths' fslext]));
        aap=aas_desc_outputs(aap, 'diffusion_session_probtrackx',[subjind sessind splitind],'lookup_tractspace_fdt_matrix2',fullfile(psesspth,'logdir',['lookup_tractspace_fdt_matrix2' fslext]));
        aap=aas_desc_outputs(aap, 'diffusion_session_probtrackx',[subjind sessind splitind],'probtrackx_log',fullfile(psesspth,'logdir','probtrackx.log' ));
        aap=aas_desc_outputs(aap, 'diffusion_session_probtrackx',[subjind sessind splitind],'waytotal',fullfile(psesspth,'logdir','waytotal' ));
            
end

