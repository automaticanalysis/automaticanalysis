function [aap, resp]=aamod_dartel_normmni(aap, task, subj)
%AAMOD_DARTEL_NORMMNISEGMENTED_MODULATED Normalise grey/white segmentations using DARTEL.
%
% After DARTEL template has been created, write out normalized
% versions (MNI space) of each subject's segmentations.
%
% This function will also do smoothing, as specified in
% aap.tasksettings.aamod_dartel_normmnisegmented.fwhm. The default
% (specified in the .xml file) is 8 mm.
%
% input streams:    dartelimported_grey
%                   dartelimported_white
%                   dartel_flowfield
%                   dartel_template
%
% output streams:   dartelnormalised_grey
%                   dartelnormalised_white


resp='';

% possible tasks 'doit','report','checkrequirements'
switch task
    case 'report'
        resp='Write smoothed normalised segmented images in MNI space for DARTEL.';
    case 'doit'
        % template
        template = aas_getfiles_bystream(aap, 'dartel_template');
        [pth,nam,ext] = fileparts(template);
        if ~exist(fullfile(aas_getsubjpath(aap,subj),[nam ext]),'file')
            copyfile(template,fullfile(aas_getsubjpath(aap,subj),[nam ext]));
        end
        template = fullfile(aas_getsubjpath(aap,subj),[nam ext]);
 
        % affine xfm
        streams=aas_getstreams(aap,'input');
        xfmi = cell_index(streams,'dartel_templatetomni_xfm');
        if xfmi && aas_stream_has_contents(aap,subj,streams{xfmi})
            xfm = load(aas_getfiles_bystream_multilevel(aap,'subject',subj,'dartel_templatetomni_xfm')); xfm = xfm.xfm;
            MMt = spm_get_space(template);
            mni.code = 'MNI152';
            mni.affine = xfm*MMt;
            [pth,nam,ext] = fileparts(template);
            save(fullfile(aas_getsubjpath(aap,subj),[nam '_2mni.mat']),'mni');            
        end
        
        % flow fields..
        job.data.subj.flowfield{1} = aas_getfiles_bystream(aap, subj, 'dartel_flowfield');
        
        % images
        imgs = '';
        streams=aap.tasklist.currenttask.outputstreams.stream;
        for streamind=1:length(streams)
            if isstruct(streams{streamind}), streams{streamind} = streams{streamind}.CONTENT; end
            if strcmp(streams{streamind},'dartel_templatetomni_xfm'), continue; end % skip
            if cell_index(aas_getstreams(aap,'input'),streams{streamind})
                imgs = strvcat(imgs, aas_getfiles_bystream(aap, subj, streams{streamind}));
            else  % renamed stream
                streamname = strrep(streams{streamind},'normalised_','');
                ind = cell_index(aas_getstreams(aap,'input'),streamname);
                if ~ind, continue; end % for dartel_templatetomni_xfm
                imgs = strvcat(imgs, aas_getfiles_bystream(aap, subj, ...
                    aap.tasklist.currenttask.inputstreams.stream{ind}));
            end
        end
        job.data.subj.images = cellstr(imgs);

        % set up job, and run
        job.template{1} = template;
        job.bb = nan(2,3);
        boundingBox = aas_getsetting(aap,'bb');
        if ~isempty(boundingBox)
            job.bb = reshape(boundingBox,2,3);
        end        
        job.vox = ones(1,3) * aap.tasklist.currenttask.settings.vox;    % voxel size
        job.fwhm = aap.tasklist.currenttask.settings.fwhm;              % smoothing
        job.preserve = aap.tasklist.currenttask.settings.preserve;      % modulation

        aas_log(aap, false, sprintf('Running with %s...', which('spm_dartel_norm_fun')));
        spm_dartel_norm_fun(job);

        % describe outputs (differ depending on modulation)
        prefix = 'w';
        if aap.tasklist.currenttask.settings.preserve, prefix = ['m' prefix]; end
        if aap.tasklist.currenttask.settings.fwhm, prefix = ['s' prefix]; end
        
        for ind=1:length(job.data.subj.images)
            [pth, nm, ext] = fileparts(job.data.subj.images{ind});
            img = fullfile(pth, [prefix nm ext]);
            aap = aas_desc_outputs(aap, subj, streams{ind}, img);
        end
        if ~exist('xfm','var')
            MMt = spm_get_space(template);
            MMm = load(fullfile(aas_getsubjpath(aap,subj),[nam '_2mni.mat']));
            MMm = MMm.mni.affine;
            xfm = MMm/MMt;
        end
        save(fullfile(aas_getsubjpath(aap,subj),'dartel_templatetomni_xfm'),'xfm')
        aap = aas_desc_outputs(aap, subj, 'dartel_templatetomni_xfm',...
            fullfile(aas_getsubjpath(aap,subj),'dartel_templatetomni_xfm.mat'));
end