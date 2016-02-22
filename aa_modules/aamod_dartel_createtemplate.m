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
        
        % retrieve external template (if any)
        template = '';
        if aas_stream_has_contents(aap,'dartel_template')
            aas_log(aap, false, sprintf('External template specified.\n'));            
            template = aas_getfiles_bystream(aap, 'dartel_template');
            aap.tasklist.currenttask.settings.exclude = 1:numel(aap.acq_details.subjects);
        end
        
        img{1} = {}; img{2} = {};
        imgTemplate{1} = {}; imgTemplate{2} = {};
        imgNoTemplate{1} = {}; imgNoTemplate{2} = {};
        for subjind = 1:numel(aap.acq_details.subjects)
            img{1}{end+1} = aas_getfiles_bystream(aap, subjind, 'dartelimported_grey');
            img{2}{end+1} = aas_getfiles_bystream(aap, subjind, 'dartelimported_white');
            if isfield(aap.tasklist.currenttask.settings,'exclude') && ...
                (...
                (isnumeric(aap.tasklist.currenttask.settings.exclude) && ...
                ~isempty(find(aap.tasklist.currenttask.settings.exclude==subjind, 1))) ||...
                (iscell(aap.tasklist.currenttask.settings.exclude) && ...
                cell_index(aap.tasklist.currenttask.settings.exclude,aas_mriname2subjname(aap.acq_details.subjects(subjind).mriname)))...
                )
                aas_log(aap,false,sprintf('Excluding subject: %s\n',aap.acq_details.subjects(subjind).mriname));
                imgNoTemplate{1}{end+1} = img{1}{end};
                imgNoTemplate{2}{end+1} = img{2}{end};
            else           
                imgTemplate{1}{end+1} = img{1}{end};
                imgTemplate{2}{end+1} = img{2}{end};
            end
        end

        % Set up job
        % Dynamic based on actual tbx_cfg_dartel
        cfg = tbx_cfg_dartel;
        cfgParam = cfg.values{2}.val{2}.val{3};
        for it = 1:numel(cfgParam.val) % for each iteration
            for par = 1:numel(cfgParam.val{1}.val) % for each parameter
                param(it).(cfgParam.val{1}.val{par}.tag) = cfgParam.val{it}.val{par}.val{1};
            end
        end
        cfgOptim = cfg.values{2}.val{2}.val{4};
        for par = 1:numel(cfgOptim.val) % for each parameter
            optim(it).(cfgOptim.val{par}.tag) = cfgOptim.val{par}.val{1};
        end
%       % Hard coded based on tbx_cfg_dartel 16 May 2012 r4667
%         param = struct(...
%             'its',{3, 3, 3, 3, 3, 3},...
%             'rparam',{[4 2 1e-6], [2 1 1e-6], [1 0.5 1e-6],...
%                       [0.5 0.25 1e-6], [0.25 0.125 1e-6], [0.25 0.125 1e-6]},...
%             'K',{0, 0, 1, 2, 4, 6},...
%             'slam',{16, 8, 4, 2, 1, 0.5});
%         optim = struct('lmreg', 0.01, 'cyc', 3, 'its', 3);

        settings = struct('template', 'Template', 'rform', aap.tasklist.currenttask.settings.rform,...
                          'param', param,...
                          'optim', optim);

        % create template
        if ~isempty(imgTemplate{1})
            spm_dartel_template(struct('images', {imgTemplate}, 'settings', settings));

            % (template in first subject)
            for t = 1:6
               template(t,:) = spm_select('fplist', fileparts(imgTemplate{1}{1}), sprintf('Template_%d',t));
            end
        end
        
        % warp excluded subjects
        if ~isempty(imgNoTemplate{1})
            settings = rmfield(settings,'template');
            for t = 1:6
               settings.param(t).template = {template(t,:)};
            end
            out = spm_dartel_warp(struct('images', {imgNoTemplate}, 'settings', settings));
            for f = 1:numel(out.files)
                movefile(out.files{f},strrep(out.files{f},'.nii','_Template.nii'),'f');
            end
        end

        % describe outputs
        aap = aas_desc_outputs(aap, 'dartel_template', template(6,:));

        % flow fields
        for subjind = 1:length(aap.acq_details.subjects)
            pth = fileparts(img{1}{subjind});
            flowimg = spm_select('fplist', pth, '^u_');
            aap = aas_desc_outputs(aap, subjind, 'dartel_flowfield', flowimg);
        end

end









