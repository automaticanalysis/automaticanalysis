% AA module - write normalised EPIs
% [aap,resp]=aamod_norm_write(aap,task,i,j)
% Rhodri Cusack MRC CBU Cambridge Jan 2006-Aug 2007
% Resamples EPIs using *_seg_sn.mat file [if present] or *_sn.mat file
% Changed domain to once per session for improved performance when parallel

function [aap,resp]=aamod_norm_write(aap,task,i,j)

resp='';

switch task
    case 'report'
    case 'doit'
        
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
            end;
            subj.imgs = strvcat(subj.imgs, P);
            % delete previous because otherwise nifti write routine doesn't
            % save disc space when you reslice to a coarser voxel
            for c=1:size(P,1)
                [pth fle ext]=fileparts(P(c,:));
                [s w]=aas_shell(['rm ' fullfile(pth,['w' fle ext])],true); % quietly
            end;
            
            % now write normalised
            if (length(subj.imgs)>0)
                spm_write_sn(subj.imgs,subj.matname,aap.spm.defaults.normalise.write);
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
