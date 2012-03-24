% AA module - write normalised mean EPI only
% Rhodri Cusack MRC CBU Cambridge Jan 2006
% Resamples EPIs using *_seg_sn.mat file [if present] or *_sn.mat file
% @@@ THIS IS NOT YET TRANSFORMED TO AA4 @@@

function [aap,resp]=aamod_norm_write_meanonly(aap,task,i)

resp='';

switch task
    case 'domain'
        resp='subject';   % this module needs to be run once per subject

    case 'description'
        resp='SPM5 write normalised (mean only)';

    case 'summary'
        subjpath=aas_getsubjpath(i);
        resp=sprintf('Write normalised mean run on %s\n',subjpath);
    case 'report'
    case 'doit'

        % get the subdirectories in the main directory
        subj_dir = aas_getsubjpath(aap,i);

        subj.imgs = [];

        % Parameter file name
        
        structdir=fullfile(aas_getsubjpath(aap,i),aap.directory_conventions.structdirname);
        % changed to s* [de 210606]
        subj.matname = dir(fullfile(structdir,['ms*' aap.acq_details.subjects_structuralfn{i} '*seg_sn.mat']));

        if (~length(subj.matname))
            subj.matname = dir(fullfile(structdir,['ms*' aap.acq_details.subjects_structuralfn{i} '*_sn.mat']));
            if (~length(subj.matname))

                aas_log(aap,1,sprintf('Cannot find normalisation _sn file in %s',fullfile(aas_getsubjpath(aap,i),aap.directory_conventions.structdirname)));
            end;
        end;

        % Add EPI mean
        P=aas_getimages(aap,i,1,['mean']);
      

        % now write normalised
        if (length(P)>0)
            spm_write_sn(P,fullfile(structdir,subj.matname.name),aap.spm.defaults.normalise.write);
        end;
    case 'checkrequirements'

    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;



