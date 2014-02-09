% AA module - write normalised EPIs useing DARTEL
% [aap,resp]=aamod_norm_write_dartel(aap,task,subj,sess)
% Rhodri Cusack MRC CBU Cambridge Jan 2006-Aug 2007
% Tibor Auer MRC CBU Cambridge 2012-2013

function [aap,resp]=aamod_norm_write_dartel(aap,task,subj,sess)

resp='';

switch task
    case 'report' % [TA]
        if ~exist('sess','var'), sess = 1; end
        % find out what streams we should normalise
		streams=aap.tasklist.currenttask.outputstreams.stream;
        for streamind=1:length(streams)
            if isstruct(streams{streamind}), streams{streamind} = streams{streamind}.CONTENT; end
            fn = ['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_' streams{streamind} '2MNI.jpg'];
            if ~exist(fullfile(aas_getsesspath(aap,subj,sess),fn),'file')
                fsl_diag(aap,subj,sess);
            end
            % Single-subjetc
            fdiag = dir(fullfile(aas_getsesspath(aap,subj,sess),'diagnostic_*.jpg'));
            for d = 1:numel(fdiag)
                aap = aas_report_add(aap,subj,'<table><tr><td>');
                aap=aas_report_addimage(aap,subj,fullfile(aas_getsesspath(aap,subj,sess),fdiag(d).name));
                aap = aas_report_add(aap,subj,'</td></tr></table>');
            end
            % Study summary
            aap = aas_report_add(aap,'reg',...
                ['Subject: ' basename(aas_getsubjpath(aap,subj)) '; Session: ' basename(aas_getsesspath(aap,subj,sess)) ]);
            aap=aas_report_addimage(aap,'reg',fullfile(aas_getsesspath(aap,subj,sess),fn));
        end
    case 'doit'
        % Is session specified in task header?
        if (isfield(aap.tasklist.currenttask.settings,'session'))
            sess = aap.tasklist.currenttask.settings.session;
        end
        
        % set up job
        % template
        template = aas_getfiles_bystream(aap, 'dartel_template');
        % flow fields..
        job.data.subj.flowfield{1} = aas_getfiles_bystream(aap, subj, 'dartel_flowfield');
        job.template{1} = template;
        job.bb = nan(2,3);
        job.vox = aap.tasklist.currenttask.settings.vox;    % voxel size
        job.fwhm = aap.tasklist.currenttask.settings.fwhm;              % smoothing
        job.preserve = aap.tasklist.currenttask.settings.preserve;      % modulation
        
		% find out what streams we should normalise
        streams=aap.tasklist.currenttask.outputstreams.stream;        
        for streamind=1:length(streams)
            imgs = [];
            % Image to reslice
            if isstruct(streams{streamind}), streams{streamind} = streams{streamind}.CONTENT; end
            if exist('sess','var')
                P = aas_getfiles_bystream(aap,subj,sess,streams{streamind});
            else
                P = aas_getfiles_bystream(aap,subj,streams{streamind});
            end
            % exclude image already normalised
            f = basename(P);
            P = P(f(:,1) ~= 'w',:);
            imgs = strvcat(imgs, P);
            % delete previous because otherwise nifti write routine doesn't
            % save disc space when you reslice to a coarser voxel
            for c=1:size(P,1)
                [pth fle ext]=fileparts(P(c,:));
                [s w] = aas_shell(['rm ' fullfile(pth,[aap.spm.defaults.normalise.write.prefix fle ext])],true); % quietly
            end;
            job.data.subj.images = cellstr(imgs);
            
            aas_log(aap, false, sprintf('Running with %s...', which('spm_dartel_norm_fun')));
            spm_dartel_norm_fun(job);
            
            % describe outputs (differ depending on modulation)
            if job.preserve==1
                prefix = 'smw';
            else
                prefix = 'sw';
            end

            % describe outputs
            wimgs=[];
            for ind=1:length(job.data.subj.images)
                [pth, nme, ext] = fileparts(job.data.subj.images{ind});
                wimgs = strvcat(wimgs,fullfile(pth,[prefix nme ext]));
            end
            if (exist('sess','var'))
                aap=aas_desc_outputs(aap,subj,sess,streams{streamind},wimgs);
            else
                aap=aas_desc_outputs(aap,subj,streams{streamind},wimgs);
            end
        end
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
end

function fsl_diag(aap,i,j) % [TA]
% Create FSL-like overview using T1 instead of the template

% find out what streams we should normalise
streams=aap.tasklist.currenttask.outputstreams.stream;
for streamind=1:length(streams)
    if isstruct(streams{streamind}), streams{streamind} = streams{streamind}.CONTENT; end
    % Obtain the structural
	sP = aas_getfiles_bystream(aap,i,'structural');
    if (strcmp(aap.tasklist.currenttask.domain,'session'))
        % Obtain the first EPI
        fP = aas_getimages_bystream(aap,i,j,streams{streamind});
    else
        fP = aas_getfiles_bystream(aap,i,streams{streamind});
    end
    [pP, fP, eP] = fileparts(fP(1,:));
    fP = fullfile(pP,['sw' fP eP]);
    % Overlays
    sess_dir = aas_getsesspath(aap,i,j);
    iP = fullfile(sess_dir,['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_' streams{streamind} '2MNI']);
    aas_runfslcommand(aap,sprintf('slices %s %s -s 2 -o %s.gif',sP,fP,iP));
    [img,map] = imread([iP '.gif']); s3 = size(img,1)/3;
    img = horzcat(img(1:s3,:,:),img(s3+1:2*s3,:,:),img(s3*2+1:end,:,:));
    imwrite(img,map,[iP '.jpg']); delete([iP '.gif']);
    iP = fullfile(sess_dir,['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_MNI2' streams{streamind}]);
    aas_runfslcommand(aap,sprintf('slices %s %s -s 2 -o %s.gif',fP,sP,iP));
    [img,map] = imread([iP '.gif']); s3 = size(img,1)/3;
    img = horzcat(img(1:s3,:,:),img(s3+1:2*s3,:,:),img(s3*2+1:end,:,:));
    imwrite(img,map,[iP '.jpg']); delete([iP '.gif']);
end
end
