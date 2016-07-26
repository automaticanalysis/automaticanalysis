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
        localpath = aas_getpath_bydomain(aap,domain,cell2mat(varargin));
        
        % Process streams
        [diagstream, mainstream] = process_streams(aap);
        d = dir(fullfile(localpath,'diagnostic_aas_checkreg_*'));
        if isempty(d)
            aas_checkreg(aap,domain,cell2mat(varargin),diagstream,'structural');
            aas_checkreg(aap,domain,cell2mat(varargin),mainstream,'structural');
        end
        subj = varargin{1};
        fdiag = dir(fullfile(localpath,'diagnostic_*.jpg'));
        for d = 1:numel(fdiag)
            aap = aas_report_add(aap,subj,'<table><tr><td>');
            imgpath = fullfile(localpath,fdiag(d).name);
            aap=aas_report_addimage(aap,subj,imgpath);
            [p f] = fileparts(imgpath); avipath = fullfile(p,[strrep(f(1:end-2),'slices','avi') '.avi']);
            if exist(avipath,'file'), aap=aas_report_addimage(aap,subj,avipath); end
            aap = aas_report_add(aap,subj,'</td></tr></table>');
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
        [diagstream, mainstream, wbstream] = process_streams(aap);
        
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
        mEPIimg = aas_getfiles_bystream(aap,domain,cell2mat(varargin),diagstream);
        if numel(spm_vol(mEPIimg)) > 1
            aas_log(aap, false, 'Found more than 1 mean functional images, using first.');
            mEPIimg = [deblank(mEPIimg(1,:)),',1'];
        end
        
        % Check local wholebrain EPI
        WBimg = '';
        if ~isempty(wbstream) && aas_stream_has_contents(aap,domain,cell2mat(varargin),wbstream)
            WBimg = aas_getfiles_bystream(aap,domain,cell2mat(varargin),wbstream);
            if size(WBimg,1) > 1
                aas_log(aap, false, sprintf('Found more than 1 wholebrain images, using %s %d', wbstream,...
                    aap.tasklist.currenttask.settings.structural));
            end
        end
        
        % Look for xfm t1totemplate
        load(aas_getfiles_bystream(aap,subj,'t1totemplate_xfm'));
        
        aas_log(aap,false,sprintf(['\tto template realignment parameters:\n' ...
            '\tx: %0.4f   y: %0.4f   z: %0.4f   p: %0.4f   r: %0.4f   j: %0.4f'], ...
            xfm(1), xfm(2), xfm(3), xfm(4), xfm(5), xfm(6)))
        
        %% 1) Mean Functional to T1 template (reorient)
        % Set the new space for the mean functional
        for img = {mEPIimg WBimg}
            if isempty(img{1}), continue; end
            spm_get_space(img{1}, spm_matrix(xfm)\spm_get_space(img{1}));
        end

        %% 2) Mean Functional to Structural (coregister)

        if ~isempty(WBimg) % two-step coreg
            % Coregister wholebrain EPI to structural
            x1 = spm_coreg(spm_vol(deblank(Simg(aap.tasklist.currenttask.settings.structural,:))), ...
                spm_vol(WBimg), ...
                flags.estimate);
            % Set the new space for the wholebrain EPI
            spm_get_space(WBimg, spm_matrix(x1)\spm_get_space(WBimg));
            
            % Coregister mean EPI to wholebrain EPI
            flags.estimate.cost_fun = 'ncc'; % whithin-modality
            x2 = spm_coreg(spm_vol(deblank(WBimg)),spm_vol(WBimg),flags.estimate);
            % Set the new space for the mean EPI
            spm_get_space(WBimg, spm_matrix(x2)\spm_get_space(WBimg));
            x = x1 + x2;
        else
            % Coregister mean EPI to structural
            x = spm_coreg(spm_vol(deblank(Simg(aap.tasklist.currenttask.settings.structural,:))), ...
                spm_vol(mEPIimg), ...
                flags.estimate);
            % Set the new space for the mean EPI
            spm_get_space(mEPIimg, spm_matrix(x)\spm_get_space(mEPIimg));
        end
                
        aas_log(aap,false,sprintf(['\tmean EPI to structural realignment parameters:\n' ...
            '\tx: %0.4f   y: %0.4f   z: %0.4f   p: %0.4f   r: %0.4f   j: %0.4f'], ...
            x(1), x(2), x(3), x(4), x(5), x(6)))
        
        %% 3) Now apply this transformation to all the EPI images
        % The mean EPI will already be in the space required for the
        % individual EPIs. Hence, we can...
        
        % Again, get space of mean functional
        MM = spm_get_space(mEPIimg(1,:));
        
        % Locate all the EPIs we want to coregister
        EPIimg = aas_getfiles_bystream(aap,domain,cell2mat(varargin),mainstream);
        for e = 1:size(EPIimg,1)
            % Apply the space of the coregistered mean EPI to the
            % remaining EPIs (safest solution!)
            spm_get_space(deblank(EPIimg(e,:)), MM);
        end
        
        %% Describe the outputs and Diagnostics
        
        if strcmp(aap.options.wheretoprocess,'localsingle')
            aas_checkreg(aap,domain,cell2mat(varargin),diagstream,'structural');
            aas_checkreg(aap,domain,cell2mat(varargin),mainstream,'structural');
        end
        
        aap = aas_desc_outputs(aap,domain,cell2mat(varargin),mainstream,EPIimg);
        
    case 'checkrequirements'
        aas_log(aap,0,'No need to trim or skull strip structural\n' );
end
end

function [diagstream, mainstream, wbstream] = process_streams(aap)
[inpstreams, inpattr] = aas_getstreams(aap,'input');
outpstreams = aas_getstreams(aap,'output');

% select non-structural diagnostic stream
diagind = cellfun(@(x) isfield(x,'diagnostic') && x.diagnostic, inpattr); diagind(cell_index(inpstreams,'structural')) = false;
if any(diagind)
    diagstream = inpstreams{diagind};
else  % auto failed --> manual  
    if cell_index(inpstreams,'meanepi') % fMRI
        diagstream = 'meanepi';
    end
    if cell_index(inpstreams,'MTI_baseline') % MTI
        diagstream = 'MTI_baseline';
    end
end

if cell_index(inpstreams,'wholebrain') % partial volume acquisition
    wbstream = inpstreams{cell_index(inpstreams,'wholebrain')};
else
    wbstream = '';
end
mainstream = outpstreams{1};
end