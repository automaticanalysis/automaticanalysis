function [aap,resp]=aamod_applyVDM(aap,task,subj,sess)

resp='';

switch task
    case 'report'
    case 'doit'
        domain = aap.tasklist.currenttask.domain;
        
        % Flags
        wrapchoice = {'ROW' 'COL' 'NYI'};
        job.roptions = aas_getsetting(aap,'roptions',sess);
        
        % EPI PE Direction
        EPI_DICOMHEADERS = load(aas_getfiles_bystream(aap,domain,[subj,sess],aas_getstreams(aap,'input',2))); EPI_DICOMHEADERS = EPI_DICOMHEADERS.DICOMHEADERS{1};
        if isfield(EPI_DICOMHEADERS,'InPlanePhaseEncodingDirection') % Siemens, GE
            job.roptions.wrap = strcmp(wrapchoice,deblank(EPI_DICOMHEADERS.InPlanePhaseEncodingDirection));
            job.roptions.pedir = find(job.roptions.wrap);
        else
            aas_log(aap,false,sprintf('WARNING: Field for phase encoding direction not found!\nWARNING: Manual setting is used: %s.',wrapchoice{job.roptions.pedir}));
            if ~isempty(job.roptions.pedir)
                job.roptions.wrap = zeros(1,3); job.roptions.wrap(job.roptions.pedir) = 1;
            else 
                aas_log(aap,true,'ERROR: No value is specified'); 
            end
        end
        
        job.data.scans = cellstr(aas_getfiles_bystream(aap,domain,[subj,sess],aas_getstreams(aap,'input',1)));
        job.data.vdmfile = cellstr(aas_getfiles_bystream(aap,domain,[subj,sess],'fieldmap'));

        FieldMap_applyvdm(job);
        
        aap=aas_desc_outputs(aap,domain,[subj,sess],aas_getstreams(aap,'output',1),spm_file(job.data.scans{1},'prefix',job.roptions.prefix));                
end
