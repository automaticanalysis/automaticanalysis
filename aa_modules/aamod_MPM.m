function [aap,resp]=aamod_MPM(aap,task,subj)

resp='';

switch task
    case 'report'
        
    case 'doit'
        
        mpmfns = cellstr(aas_getfiles_bystream(aap,'special_session',[subj,cell_index({aap.acq_details.special_sessions.name},'MPM')],'MPM'));
        
        job.b1_type = '3D_EPI_v2b';
        
        job.data_spec.subj.raw_fld.b0 = [...
            mpmfns(cell_index(mpmfns,'serie01')),...
            {''},...
            mpmfns(cell_index(mpmfns,'serie02')),...
            ];
        job.data_spec.subj.raw_fld.b1 = mpmfns(cell_index(mpmfns,'serie03'))';
        job.data_spec.subj.raw_mpm.MT = mpmfns(cell_index(mpmfns,'serie04'))';
        job.data_spec.subj.raw_mpm.PD = mpmfns(cell_index(mpmfns,'serie05'))';
        job.data_spec.subj.raw_mpm.T1 = mpmfns(cell_index(mpmfns,'serie06'))';
        
        aas_makedir(aap,fullfile(aas_getsubjpath(aap,subj),'MPM'));
        job.data_spec.subj.output.outdir = cellstr(fullfile(aas_getsubjpath(aap,subj),'MPM'));
        
        out = vbq_mpr_b0_b1(job);
        
        aap=aas_desc_outputs(aap,'subject',subj,'mt',out.MT{1});
        aap=aas_desc_outputs(aap,'subject',subj,'mtr',spm_select('FPListRec',spm_file(job.data_spec.subj.raw_mpm.MT{1},'path'),'.*_MTR.nii$'));
        aap=aas_desc_outputs(aap,'subject',subj,'mtw',spm_select('FPListRec',spm_file(job.data_spec.subj.raw_mpm.MT{1},'path'),'.*_MTw.nii$'));
        aap=aas_desc_outputs(aap,'subject',subj,'pd',strrep(out.A{1},'_A.nii','_PDw.nii'));
        aap=aas_desc_outputs(aap,'subject',subj,'r1',out.R1{1});
        aap=aas_desc_outputs(aap,'subject',subj,'r2star',out.R2s{1});
        aap=aas_desc_outputs(aap,'subject',subj,'r2star_ols',spm_select('FPListRec',spm_file(job.data_spec.subj.raw_mpm.MT{1},'path'),'.*_R2s_OLS.nii$'));
        aap=aas_desc_outputs(aap,'subject',subj,'t1',out.T1w{1});
        
    case 'checkrequirements'
        
end
