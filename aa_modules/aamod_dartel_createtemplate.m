function [aap,resp] = aamod_dartel_createtemplate(aap,task)
%AAMOD_DARTEL_CREATETEMPLATE Create a template using DARTEL.
%
% Template creation can be quite time consuming; for example, 48
% hours for a template of ~80 subjects is not unreasonable.
%
% input streams:    dartelimported_grey
%                   dartelimported_white
%
% output streams:   dartel_template
%                   dartel_flowfield

resp='';

switch task
    case 'domain'
        resp='study';
    case 'description'
        resp='Create DARTEL template';
    case 'doit'
        % report which version of spm_dartel_template we are using
        aas_log(aap, false, sprintf('Using %s.\n', which('spm_dartel_template')));

        % number of tissue classes included
        numtissues = 2;

        for subjind = 1:length(aap.acq_details.subjects)
            img{1}{subjind} = aas_getfiles_bystream(aap, subjind, 'dartelimported_grey');
            img{2}{subjind} = aas_getfiles_bystream(aap, subjind, 'dartelimported_white');
        end


        % Set up job
        % below based on tbx_cfg_dartel 16 May 2012 r4667
        param = struct(...
            'its',{3, 3, 3, 3, 3, 3},...
            'rparam',{[4 2 1e-6], [2 1 1e-6], [1 0.5 1e-6],...
                      [0.5 0.25 1e-6], [0.25 0.125 1e-6], [0.25 0.125 1e-6]},...
            'K',{0, 0, 1, 2, 4, 6},...
            'slam',{16, 8, 4, 2, 1, 0.5});

        settings = struct('template', 'Template', 'rform', aap.tasklist.currenttask.settings.rform,...
                          'param', param,...
                          'optim', struct('lmreg', 0.01, 'cyc', 3, 'its', 3));

        % create the job
        job = struct('images', {img}, 'settings', settings);

        % run the script
        spm_dartel_template(job);

        % describe outputs

        % (template in first subject)
        [pth, nm] = fileparts(img{1}{1});
        templateimg = spm_select('fplist', pth, 'Template_6');
        aap = aas_desc_outputs(aap, 'dartel_template', templateimg);

        % flow fields
        for subjind = 1:length(aap.acq_details.subjects)
            [pth, nm] = fileparts(img{1}{subjind});
            flowimg = spm_select('fplist', pth, '^u_');
            aap = aas_desc_outputs(aap, subjind, 'dartel_flowfield', flowimg);
        end

end









