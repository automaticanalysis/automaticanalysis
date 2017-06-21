function [aap,resp]=aamod_stringentwbmask_secondlevel(aap,task,i,j)
% AA module - create a mean structural for your sample. 
% 19/2/2010 J Carlin

resp='';

switch task
    case 'doit'

		% For each subject
		for sub = 1:length(aap.acq_details.subjects)
			% Get T1
            p = aas_getfiles_bystream(aap,sub,'structural');
            % for some reason the true normalised image is now the second
            % here...
            V(sub) = spm_vol(p(2,:));
		end

		% Load all those vols
		xyz = spm_read_vols(V);

		V_out = V(end);
        outdir = fullfile(aas_getstudypath(aap),'structurals');
        mkdirifneeded(outdir);
		V_out.fname = fullfile(outdir,'mean_structural.nii');
		spm_write_vol(V_out,nanmean(xyz,4));
        aap = aas_desc_outputs(aap,'structural_group',V_out.fname);

        % for diagnostic purposes, also record variance - useful for
        % identifying poorly normalised or variable regions
        V_out.fname = fullfile(outdir,'std_structural.nii');
        spm_write_vol(V_out,nanstd(xyz,[],4));

    case 'checkrequirements'

    otherwise
     aas_log(aap,1,sprintf('Unknown task %s',task));
end



