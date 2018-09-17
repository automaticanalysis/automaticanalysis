% AA module - runs FSL POSSUM
% i.e. Simulated fMRI images!
% Rhodri Cusack UWO Jun 2011
% rhodri@cusacklab.org

function [aap,resp]=aamod_possum(aap,task,job)

resp='';

switch task
    case 'report'
        
    case 'getparallelparts'
        resp=[];
        for pp=1:20
            resp{pp}.id=sprintf('part%d',pp);
            resp{pp}.timeinseconds=aap.tasklist.currenttask.settings.pulse_commandlineparameters.tr*(pp-1);
            resp{pp}.ppnum=pp;
        end;
    case 'doit'
        aas_log(aap,0,sprintf('Now trying to deal with %s',job.id));
        
        outputpath=fullfile(aas_getstudypath(aap),job.id);
        aas_makedir(aap,outputpath);
        aas_makedir(aap,fullfile(outputpath,'diff_proc'));
        
        % Set up command for fsl tool "pulse"
        cmd=['pulse -i ' aap.tasklist.currenttask.settings.inputbrain ' -o ' outputpath '/pulse -v '];
        clp=aap.tasklist.currenttask.settings.pulse_commandlineparameters;
        clp_fld=fieldnames(clp);
        for fldind=1:length(clp_fld)
            val=clp.(clp_fld{fldind});
            if (isnumeric(val))
                val=num2str(val);
            end;
            cmd=[cmd ' --' clp_fld{fldind} '=' val];
        end;
        
        % Run pulse
        [s w]=aas_runfslcommand(aap,cmd);
        
        % Set up command for fsl tool "possum"
        %            {aap.tasklist.currenttask.settings.motion,'motion'}, ...
        
        filestocopy={{aap.tasklist.currenttask.settings.inputbrain,'brain.nii.gz'}, ...
            {aap.tasklist.currenttask.settings.sliceprofile,'slcprof'}, ...
            {aap.tasklist.currenttask.settings.motion,'motion_canned'}, ...
            {aap.tasklist.currenttask.settings.mrparameters,'MRpar'}};
        for fnind=1:length(filestocopy)
            % Use cp command with aas_runfslcommand, so that ${FSLDIR} in
            % paths will be expanded as needed
            cmd=['cp ' filestocopy{fnind}{1} ' ' fullfile(outputpath,filestocopy{fnind}{2})];
            [s w]=aas_runfslcommand(aap,cmd);
        end;
        
        
        
        % Repeatable random sequence
        RandStream.setDefaultStream(RandStream('mt19937ar','seed',1234))
        
        % Position and rotations are random walk
        
        nsteps=200;
        posstep=[0.001 0.001 0.001  0.01 0.01 0.01]/sqrt(nsteps);
        pos=zeros(1,6);
        
        % Recreate random walk for previous volumes to get motions for current volume
        for pp=1:job.ppnum
            motionlist=[];
            for time=linspace(0,0.95*aap.tasklist.currenttask.settings.pulse_commandlineparameters.tr,nsteps);
                motionlist=[motionlist;time pos];
                pos=pos+(rand(1,6)*2-1).*posstep;
            end;
        end;
        
        % Dump motion list
        fid=fopen(fullfile(outputpath,'motion'),'w');
        for rowind=1:size(motionlist,1)
            fprintf(fid,'   %f',[motionlist(rowind,:)]);
            fprintf(fid,'\n');
        end;
        fclose(fid);
        
        % Now run possum
        cmd=['possumX ' outputpath]
        [s w]=aas_runfslcommand(aap,cmd);
        
        
        aap=aas_desc_outputs(aap,sprintf('image_abs_%s',job.id),fullfile(outputpath,'image_abs.nii.gz'));
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
end
