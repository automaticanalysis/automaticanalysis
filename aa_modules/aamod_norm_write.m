% AA module - write normalised EPIs
% [aap,resp]=aamod_norm_write(aap,task,i,j)
% Rhodri Cusack MRC CBU Cambridge Jan 2006-Aug 2007
% Resamples EPIs using *_seg_sn.mat file [if present] or *_sn.mat file
% Changed domain to once per session for improved performance when parallel
% Tibor Auer MRC CBU Cambridge 2012-2013

function [aap,resp]=aamod_norm_write(aap,task,i,j)

resp='';

switch task
	case 'report' % [TA]
        if ~exist('j','var'), j = 1; end
        fn = ['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_epi2MNI.jpg'];
        if ~exist(fullfile(aas_getsesspath(aap,i,j),fn),'file')
            fsl_diag(aap,i,j);
        end
        % Single-subjetc
        fdiag = dir(fullfile(aas_getsesspath(aap,i,j),'diagnostic_*.jpg'));
        for d = 1:numel(fdiag)
            aap = aas_report_add(aap,i,'<table><tr><td>');
            aap=aas_report_addimage(aap,i,fullfile(aas_getsesspath(aap,i,j),fdiag(d).name));
            aap = aas_report_add(aap,i,'</td></tr></table>');
        end
        % Study summary
        aap = aas_report_add(aap,'reg',...
            ['Subject: ' basename(aas_getsubjpath(aap,i)) '; Session: ' basename(aas_getsesspath(aap,i,j)) ]);
        aap=aas_report_addimage(aap,'reg',fullfile(aas_getsesspath(aap,i,j),fn));
    case 'doit'
        
        voxelSize = aap.tasklist.currenttask.settings.vox; % in case we want something other than default voxel size
        
        % get the subdirectories in the main directory
        subj_dir = aas_getsubjpath(aap,i); 
        
        % get sn mat file from normalisation
        subj.matname = aas_getfiles_bystream(aap,i,'normalisation_seg_sn');
        
        streams=aap.tasklist.currenttask.inputstreams;
        
        % find out what streams we should normalise
        streams=streams.stream(~[strcmp('normalisation_seg_sn',streams.stream)]);
        
        % Is session specified in task header (used for meanepi, which only
        % occurs in session 1
        if (isfield(aap.tasklist.currenttask.settings,'session'))
            j=aap.tasklist.currenttask.settings.session;
        end;
        
        for streamind=1:length(streams)
            subj.imgs = [];
                        
            % Image to reslice
            if (exist('j','var'))
                P = aas_getfiles_bystream(aap,i,j,streams{streamind});
            else
                P = aas_getfiles_bystream(aap,i,streams{streamind});
            end
            
            
            subj.imgs = strvcat(subj.imgs, P);
            % delete previous because otherwise nifti write routine doesn't
            % save disc space when you reslice to a coarser voxel
            for c=1:size(P,1)
                [pth fle ext]=fileparts(P(c,:));
                [s w]=aas_shell(['rm ' fullfile(pth,['w' fle ext])],true); % quietly
            end;
            
            
            % set defaults
            flags = aap.spm.defaults.normalise.write;
            flags.vox = voxelSize; 
            
            % now write normalised
            if (length(subj.imgs)>0)
                spm_write_sn(subj.imgs,subj.matname, flags);
            end;
            wimgs=[];
            
            % describe outputs
            for fileind=1:size(subj.imgs,1)
                [pth nme ext]=fileparts(subj.imgs(fileind,:));
                wimgs=strvcat(wimgs,fullfile(pth,['w' nme ext]));
            end;
            if (exist('j','var'))
                aap=aas_desc_outputs(aap,i,j,streams{streamind},wimgs);
            else
                aap=aas_desc_outputs(aap,i,streams{streamind},wimgs);
            end;
        end;
        
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
end

function fsl_diag(aap,i,j) % [TA]
% Create FSL-like overview using T1 instead of the template

% Obtain structural from aamod_norm_noss
norm_dir = aas_getsubjpath(aap,i,cell_index({aap.tasklist.main.module.name},'aamod_norm_noss'));
structdir=fullfile(norm_dir,aap.directory_conventions.structdirname);
sP = dir( fullfile(structdir,[aap.spm.defaults.normalise.write.prefix 's' aap.acq_details.subjects(i).structuralfn '*.nii']));
sP = fullfile(structdir,sP(1).name);

% Obtain the first EPI
sess_dir = aas_getsesspath(aap,i,j);
fP = aas_getimages_bystream(aap,i,j,aap.tasklist.currenttask.inputstreams.stream{2});
[pP fP eP] = fileparts(fP(1,:));
fP = fullfile(pP,[aap.spm.defaults.normalise.write.prefix fP eP]);

% Overlays
iP = fullfile(sess_dir,['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_epi2MNI']);
system(sprintf('slices %s %s -s 2 -o %s.gif',sP,fP,iP));
[img,map] = imread([iP '.gif']); s3 = size(img,1)/3;
img = horzcat(img(1:s3,:,:),img(s3+1:2*s3,:,:),img(s3*2+1:end,:,:));
imwrite(img,map,[iP '.jpg']); delete([iP '.gif']);
iP = fullfile(sess_dir,['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_MNI2epi']);
system(sprintf('slices %s %s -s 2 -o %s.gif',fP,sP,iP));
[img,map] = imread([iP '.gif']); s3 = size(img,1)/3;
img = horzcat(img(1:s3,:,:),img(s3+1:2*s3,:,:),img(s3*2+1:end,:,:));
imwrite(img,map,[iP '.jpg']); delete([iP '.gif']);
end
