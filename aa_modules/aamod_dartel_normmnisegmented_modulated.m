function [aap,resp]=aamod_dartel_normmnisegmented_modulated(aap, task)
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
    case 'domain'
        resp='subject';
    case 'report'
        resp='Write smoothed normalised segmented images in MNI space for DARTEL.'
    case 'doit'
        % template
        template = aas_getfiles_bystream(aap, 'dartel_template');

        % initialize images
        allimages = {};

        for subjind = 1:length(aap.acq_details.subjects)

            % flow fields..
            job.data.subj(subjind).flowfield{1} = aas_getfiles_bystream(aap, subjind, 'dartel_flowfield');

            % images
            imgs1 = aas_getfiles_bystream(aap, subjind, 'dartelimported_grey');
            imgs2 = aas_getfiles_bystream(aap, subjind, 'dartelimported_white');

            job.data.subj(subjind).images = cellstr(strvcat(imgs1, imgs2));
        end % going through subjects

        % set up job
        job.template{1} = template;
        job.bb = nan(2,3);
        job.vox = ones(1,3) * aap.tasklist.currenttask.settings.vox;
        job.fwhm = aap.tasklist.currenttask.settings.fwhm;
        job.preserve = 1; % modulated

        aas_log(aap, false, sprintf('Running with %s...', which('spm_dartel_norm_fun')));
        spm_dartel_norm_fun(job);

        % describe outputs
        for subjind = 1:length(aap.acq_details.subjects)
            [pth, nm, ext] = fileparts(job.data.subj(subjind).images{1});
            greyimg = fullfile(pth, ['smw' nm ext]);
            aap = aas_desc_outputs(aap, subjind, 'dartelnormalised_grey', greyimg);

            [pth, nm, ext] = fileparts(job.data.subj(subjind).images{2});
            whiteimg = fullfile(pth, ['smw' nm ext]);
            aap = aas_desc_outputs(aap, subjind, 'dartelnormalised_white', whiteimg);
        end
end
