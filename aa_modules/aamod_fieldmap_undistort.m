% AA module - realignment and unwarp
% Realignment using SPM5
% i=subject num
% Rhodri Cusack MRC CBU 2004-6

function [aap,resp]=aamod_realign(aap,task,i)

resp='';

switch task
    case 'domain'
        resp='subject';
    case 'description'
        resp='Fieldmap undistort v4';
    case 'summary'
        resp='Done Fieldmap undistort v4\n';
    case 'report'

    case 'doit'
        % Get mean
        epimeanfilter=['mean'];
        epimean= aas_getimages(aap,i,1,epimeanfilter);
        % Get all EPIs to undistort
        epistoundistort=[];
        for j = aap.acq_details.selected_sessions
            tmp= aas_getimages(aap,i,j,aap.tasklist.currenttask.epiprefix);
            epistoundistort=strvcat(epistoundistort,tmp);
        end
        % Get phase and magnitude
        subjpath=aas_getsubjpath(aap,i);
        fieldmapsdir=fullfile(subjpath,aap.directory_conventions.fieldmapsdirname);
        magdir=fullfile(fieldmapsdir,'rawmag');
        magfn=dir(fullfile(magdir,'*nii'));
        if (length(magfn)<1)
            aas_log(aap,1,sprintf('No fieldmap magnitude found in %s',magdir));
        end;
        phasedir=fullfile(fieldmapsdir,'rawphase');
        phasefn=dir(fullfile(phasedir,'*nii'));
        if (length(phasefn)<1)
            aas_log(aap,1,sprintf('No fieldmap phase found in %s',phasedir));
        end;
        magfn=fullfile(magdir,magfn(1).name);
        phasefn=fullfile(phasedir,phasefn(1).name);

        % Do undistortion
        if (aap.options.fieldmapundistortversion=='fieldmap_undistort_v401')
            feval(aap.options.fieldmapundistortversion,magfn,phasefn,epimean,epistoundistort,fieldmapsdir);
        else
            subj_dir = aas_getsubjpath(aap,i);
            % Get structural for this subject
            structdir=fullfile(subj_dir,aap.directory_conventions.structdirname);
            structfn = dir( fullfile(structdir,['s' aap.acq_details.subjects_structuralfn{i} '*.nii']));
            structfn = fullfile(structdir,structfn(1).name);
            feval(aap.options.fieldmapundistortversion,magfn,phasefn,epimean,epistoundistort,fieldmapsdir,structfn);
        end;
        % Save graphical output
        figure(spm_figure('FindWin', 'Graphics'));
        print('-djpeg','-r75',fullfile(aas_getsubjpath(aap,i),'diagnostic_aamod_fieldmapundistort'));

    case 'checkrequirements'

    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
