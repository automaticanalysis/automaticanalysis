% AA module - coregister structural to mean EPI
% Coregistration of structural to mean EPI output by realignment
% Does not require skull stripping any more
% Modified for sparse imaging since prefix for mean is different
% i=subject num
% Rhodri Cusack MRC CBU 2004-6 based on original by Matthew Brett
% 
% Major changes Aug 2010: removed support for central store of structrual
% images. This code was very long in tooth, and unloved.

function [aap,resp]=aamod_coreg(aap,task,i)

resp='';

switch task
    case 'doit'
        global defaults;
        flags = defaults.coreg;
        % check local structural directory exists
        subjpath=aas_getsubjpath(aap,i);
        structdir=fullfile(subjpath,aap.directory_conventions.structdirname);
        if (~length(dir(structdir)))
            [s w]=aas_shell(['mkdir ' structdir]);
            if (s)
                aas_log(aap,1,sprintf('Problem making directory%s',structdir));
            end;
        end;
        
        % dirnames,
        % get the subdirectories in the main directory
        dirn = aas_getsesspath(aap,i,1);
        % get mean EPI stream
        PG = aas_getimages_bystream(aap,i,1,'meanepi');
        VG = spm_vol(PG);
        
        % Get path to structural for this subject
        inStream = aap.tasklist.currenttask.inputstreams.stream{1};
        structImg = aas_getfiles_bystream(aap, subjInd, inStream);                
        VF = spm_vol(structImg);

        % do coregistration
        x  = spm_coreg(VG, VF,flags.estimate);
        
        M  = inv(spm_matrix(x));
          
        spm_get_space(structfn, M*spm_get_space(structfn));
       
        aap = aas_desc_outputs(aap,i,'structural',structfn);

        % Save graphical output - this will now be done by report task
        try figure(spm_figure('FindWin', 'Graphics')); catch; figure(1); end;            
        print('-djpeg','-r75',fullfile(aas_getsubjpath(aap,i),'diagnostic_aamod_coreg'));
        
    case 'checkrequirements'
        aas_log(aap,0,'No need to trim or skull strip structural\n' );
end
