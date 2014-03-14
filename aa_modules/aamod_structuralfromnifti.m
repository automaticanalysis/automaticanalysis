% AA module - structural from NIFTI

function [aap,resp]=aamod_structuralfromnifti(aap,task,subj)

resp='';

switch task
    case 'report'
    case 'doit'
        if ~iscell(aap.acq_details.subjects(subj).structural)
            aas_log(aap,true,'Was exepcting list of filenames in cell array');
        end;
        niftifile = '';
        niftisearchpth=aas_findvol(aap,'');
        if ~isempty(niftisearchpth)
            niftifile = fullfile(niftisearchpth,aap.acq_details.subjects(subj).structural{1});
        end
        V=spm_vol(niftifile);
        sesspth=fullfile(aas_getsubjpath(aap,subj),aap.directory_conventions.structdirname);
        aas_makedir(aap,sesspth);
        [pth fle ext]=fileparts(aap.acq_details.subjects(subj).structural{1});
        Y=spm_read_vols(V(1));
        fn=fullfile(sesspth,[sprintf('%s_%04d',fle,1) '.nii']);
        V(1).fname=fn;
        V(1).n=[1 1];
        % Write out the files, now likely in 3d
        spm_write_vol(V(1),Y);        
        aap=aas_desc_outputs(aap,subj,'structural',fn);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end