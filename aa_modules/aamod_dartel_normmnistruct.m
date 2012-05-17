function [aap, resp] = aamod_dartel_normmnistruct(aap, task)
%AAMOD_DARTEL_NORMMNISTRUCT Write DARTEL-normalised structural image.
%
% After DARTEL template has been created, write out normalized
% versions (MNI space) of each subject's structural image.
%
% This function will also do smoothing, as specified in
% aap.tasksettings.aamod_dartel_normmnistruct.fwhm. The default
% (specified in the .xml file) is .5 mm (to help with aliasing).
%
% input streams:    dartel_template
%                   dartel_flowfield
%                   structural
%
% output streams:   dartelnormalised_structural

resp='';

% possible tasks 'doit','report','checkrequirements'
switch task
    case 'domain'
        resp='subject';
    case 'report'
        resp='Write smoothed normalised structural images in MNI space for DARTEL.'
    case 'doit'

        % template
        template = aas_getfiles_bystream(aap, 'dartel_template');

        % initialize images
        allimages = {};

        for subjind = 1:length(aap.acq_details.subjects)

            % flow fields..
            job.data.subj(subjind).flowfield{1} = aas_getfiles_bystream(aap, subjind, 'dartel_flowfield');

            % images
            img = aas_getfiles_bystream(aap, subjind, 'structural');

            job.data.subj(subjind).images = cellstr(img);
        end % going through subjects

        % set up job
        job.template{1} = template;
        job.bb = nan(2,3);
        job.vox = ones(1,3) * aap.tasklist.currenttask.settings.vox;
        job.fwhm = aap.tasklist.currenttask.settings.fwhm;
        job.preserve = aap.tasklist.currenttask.settings.preserve; % unmodulated/modulated

        aas_log(aap, false, sprintf('Running with %s...', which('spm_dartel_norm_fun')));
        spm_dartel_norm_fun(job);

        % describe outputs
        for subjind = 1:length(aap.acq_details.subjects)
            [pth, nm] = fileparts(job.data.subj(subjind).images{1});
            structimg = spm_select('fplist', pth, ['^smw' nm]);
            aap = aas_desc_outputs(aap, subjind, 'dartelnormalised_structural', structimg);
        end
end