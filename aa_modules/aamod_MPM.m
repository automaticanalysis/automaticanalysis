function [aap,resp]=aamod_MPM(aap,task,subj)

resp='';

switch task
    case 'report'
        
    case 'doit'
        
        mpmfns = cellstr(aas_getfiles_bystream(aap,'special_session',[subj,cell_index({aap.acq_details.special_sessions.name},'MPM')],'MPM'));
        b0fn = cellstr(aas_getfiles_bystream(aap,'special_session',[subj,cell_index({aap.acq_details.special_sessions.name},'MPM')],'fieldmap'));
        
        switch aas_getsetting(aap,'sensitivity')
            case 'RF_none'
                job.subj.sensitivity.RF_none = 'noRF';
        end
        switch aas_getsetting(aap,'b1_type')
            case 'i3D_EPI'
                job.subj.b1_type.i3D_EPI.b1input = mpmfns(cell_index(mpmfns,'serie01'));
                job.subj.b1_type.i3D_EPI.b0input = b0fn([1 2 2]);
                job.subj.b1_type.i3D_EPI.b1parameters.b1metadata = 'yes';
        end
        
        job.subj.raw_mpm.MT = mpmfns(cell_index(mpmfns,'serie02'))';
        job.subj.raw_mpm.PD = mpmfns(cell_index(mpmfns,'serie03'))';
        job.subj.raw_mpm.T1 = mpmfns(cell_index(mpmfns,'serie04'))';
        
        job.subj.output.outdir = cellstr(fullfile(aas_getsubjpath(aap,subj),'MPM'));
        
        out = hmri_run_create(job);
        
        aap=aas_desc_outputs(aap,'subject',subj,'mt',out.MT{1});
        aap=aas_desc_outputs(aap,'subject',subj,'mtw',out.MTw{1});
        aap=aas_desc_outputs(aap,'subject',subj,'pd',out.A{1});
        aap=aas_desc_outputs(aap,'subject',subj,'pdw',out.PDw{1});
        aap=aas_desc_outputs(aap,'subject',subj,'r1',out.R1{1});
        aap=aas_desc_outputs(aap,'subject',subj,'r2star',out.R2s{1});
        aap=aas_desc_outputs(aap,'subject',subj,'t1',out.T1w{1});
        
    case 'checkrequirements'
        if isempty(which('hmri_run_create'))
            aas_log(aap,true,sprintf('ERROR: hMRI toolbox not found!\n\tMake sure it is added to aap.directory_conventions.spmtoolsdir!'));
        end
end
