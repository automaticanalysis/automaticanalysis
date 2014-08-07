% AA module - extended coregistration
% Coregistration of structural to T1 template
% 1) Reorient MeanEPI to structural(output 'structural' from aamod_coreg_extended_1)
% 2) Coregister MeanEPI to structural (using xfm mat called
% 't1totemplate_xfm from aamod_coreg_extended_1)
% 3) Apply all (1 and 2) to epi.

function [aap,resp]=aamod_coreg_extended_2(aap,task,varargin)

resp='';

switch task
    case 'report' % [TA]
        domain = aap.tasklist.currenttask.domain;
        localpath = aas_getpath_bydomain(aap,domain,varargin{:});
        
        % Process streams
        diagstream = process_streams(aap);
        d = dir(fullfile(localpath,['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_structural2*']));
        if isempty(d)
            aas_checkreg(aap,domain,varargin{:},diagstream,'structural');
        end
        if strcmp(domain,'subject')
            subj = varargin{1};
            fdiag = dir(fullfile(localpath,'diagnostic_*.jpg'));
            for d = 1:numel(fdiag)
                aap = aas_report_add(aap,subj,'<table><tr><td>');
                aap=aas_report_addimage(aap,subj,fullfile(localpath,fdiag(d).name));
                aap = aas_report_add(aap,subj,'</td></tr></table>');
            end
        end
    case 'doit'
        
        global defaults
        flags = defaults.coreg;
        if isfield(aap.tasklist.currenttask.settings,'eoptions')
            fields = fieldnames(aap.tasklist.currenttask.settings.eoptions);
            for f = 1:numel(fields)
                if ~isempty(aap.tasklist.currenttask.settings.eoptions.(fields{f}))
                    flags.estimate.(fields{f}) = aap.tasklist.currenttask.settings.eoptions.(fields{f});
                end
            end
        end
        
        %% 0) Check that the templates and images we need exist!
        subj = varargin{1};
        
        % Process streams
        domain = aap.tasklist.currenttask.domain;
        [diagstream, mainstream] = process_streams(aap);
        
        % Get the T1 template
        sTimg = fullfile(spm('dir'), aap.directory_conventions.T1template);
        if ~exist(sTimg, 'file')
            aas_log(aap, true, sprintf('Couldn''t find template T1 image %s.', sTimg));
        end
        
        % Check local structural
        Simg = aas_getfiles_bystream(aap,subj,'structural');
        if size(Simg,1) > 1
            aas_log(aap, false, sprintf('Found more than 1 structural images, using structural %d', ...
                aap.tasklist.currenttask.settings.structural));
        end
        
        % Look for mean functional
        mEPIimg = aas_getfiles_bystream(aap,domain,varargin{:},diagstream);
        if size(mEPIimg,1) > 1
            aas_log(aap, false, 'Found more than 1 mean functional images, using first.');
            mEPIimg = deblank(mEPIimg(1,:));
        end
        
        % Look for xfm t1totemplate
        load(aas_getfiles_bystream(aap,subj,'t1totemplate_xfm'));
        Me = inv(spm_matrix(xfm));
        
        fprintf(['\tto template realignment parameters:\n' ...
            '\tx: %0.4f   y: %0.4f   z: %0.4f   p: %0.4f   r: %0.4f   j: %0.4f\n'], ...
            xfm(1), xfm(2), xfm(3), xfm(4), xfm(5), xfm(6))
        
        %% 1) Mean Functional to T1 template (reorient)
        % Set the new space for the mean functional
        spm_get_space(mEPIimg, Me*spm_get_space(mEPIimg));

        %% 2) Mean Functional to Structural (coregister)
        
        % Coregister mean EPI to structural
        x = spm_coreg(spm_vol(deblank(Simg(aap.tasklist.currenttask.settings.structural,:))), ...
            spm_vol(mEPIimg), ...
            flags.estimate);
        % Set the new space for the mean EPI
        spm_get_space(mEPIimg, spm_matrix(x)\spm_get_space(mEPIimg));
                
        fprintf(['\tmean EPI to structural realignment parameters:\n' ...
            '\tx: %0.4f   y: %0.4f   z: %0.4f   p: %0.4f   r: %0.4f   j: %0.4f\n'], ...
            x(1), x(2), x(3), x(4), x(5), x(6))
        
        %% 3) Now apply this transformation to all the EPI images
        % The mean EPI will already be in the space required for the
        % individual EPIs. Hence, we can...
        
        % Again, get space of mean functional
        MM = spm_get_space(mEPIimg(1,:));
        
        % Locate all the EPIs we want to coregister
        EPIimg = aas_getfiles_bystream(aap,domain,varargin{:},mainstream);
        for e = 1:size(EPIimg,1)
            % Apply the space of the coregistered mean EPI to the
            % remaining EPIs (safest solution!)
            spm_get_space(deblank(EPIimg(e,:)), MM);
        end
        
        %% Describe the outputs and Diagnostics
        
        aap = aas_desc_outputs(aap,domain,varargin{:},diagstream,mEPIimg);
        if strcmp(aap.options.wheretoprocess,'localsingle')
            aas_checkreg(aap,domain,varargin{:},diagstream,'structural');
        end
        
        aap = aas_desc_outputs(aap,domain,varargin{:},mainstream,EPIimg);
        
    case 'checkrequirements'
        aas_log(aap,0,'No need to trim or skull strip structural\n' );
end
end

function [diagstream, mainstream] = process_streams(aap)
outpstreams = aas_getstreams(aap,'out');
if cell_index(outpstreams,'rois') % fMRI
    diagstream = 'meanepi';
    mainstream = 'rois';
end
end