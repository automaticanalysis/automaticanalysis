% General coregistration
% Target: first (existing) inputstream
% Source: second
% Other (optional): third, etc.

function [aap,resp]=aamod_coreg_general(aap,task,varargin)

resp='';

switch task
    case 'report'
        domain = aap.tasklist.currenttask.domain;
        localpath = aas_getpath_bydomain(aap,domain,varargin{:});
        
        % Process streams
        inpstream = aap.tasklist.currenttask.settings.inputstreams.stream;
        outpstream = aap.tasklist.currenttask.settings.outputstreams.stream;
        if ~iscell(outpstream), outpstream = cellstr(outpstream); end;
        tInd = 1;
        while ~aas_stream_has_contents(aap,aap.tasklist.currenttask.domain,varargin{:},inpstream{tInd})
            tInd = tInd + 1;
        end
        targetimfn = aas_getfiles_bystream(aap,domain,varargin{:},inpstream{tInd}); % occasional "." in streamname may cause problem for aas_checkreg
        
        d = dir(fullfile(localpath,['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_*']));
        if isempty(d)
            aas_checkreg(aap,domain,varargin{:},inpstream{tInd+1},targetimfn);
            if numel(inpstream) > (tInd+1)
                for s = 1:numel(inpstream)-(tInd+1)
                    aas_checkreg(aap,domain,varargin{:},inpstream{tInd+1+s},targetimfn);
                end
            end
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
        
        %% Init
        flags = aap.spm.defaults.coreg;
        % update flags
        flags.write.which = [1 0]; % do not (re)write target and mean
        if isfield(aap.tasklist.currenttask.settings,'eoptions')
            fields = fieldnames(aap.tasklist.currenttask.settings.eoptions);
            for f = 1:numel(fields)
                if ~isempty(aap.tasklist.currenttask.settings.eoptions.(fields{f}))
                    flags.estimate.(fields{f}) = aap.tasklist.currenttask.settings.eoptions.(fields{f});
                end
            end
        end
        if isfield(aap.tasklist.currenttask.settings,'roptions')
            fields = fieldnames(aap.tasklist.currenttask.settings.eoptions);
            for f = 1:numel(fields)
                if ~isempty(aap.tasklist.currenttask.settings.eoptions.(fields{f}))
                    flags.write.(fields{f}) = aap.tasklist.currenttask.settings.eoptions.(fields{f});
                end
            end
        end

        subj = varargin{1};
        domain = aap.tasklist.currenttask.domain;
        
        %% Data
        % Process streams
        inpstream = aap.tasklist.currenttask.settings.inputstreams.stream;
        outpstream = aap.tasklist.currenttask.settings.outputstreams.stream;
        if ~iscell(outpstream), outpstream = cellstr(outpstream); end;
        
        % Get target image:        
        tInd = 1;
        while ~aas_stream_has_contents(aap,aap.tasklist.currenttask.domain,varargin{:},inpstream{tInd})
            tInd = tInd + 1;
        end
        fprintf('\nTarget stream: %s\n',inpstream{tInd});
        targetimfn = aas_getfiles_bystream(aap,domain,varargin{:},inpstream{tInd});
        if size(targetimfn,1) > 1 % multiple images (possibly after normalisation)
            fns = spm_file(targetimfn,'basename');
            targetimfn = targetimfn(fns(:,1)~='w',:); % omit normalised
        end
        if size(targetimfn,1) > 1 % multiple images --> ?
            aas_log(aap,true,sprintf('ERROR: multiple target image found:%s',targetimfn'));
        end
        
        % Get image to coregister ('source'):
        sourceimfn = aas_getfiles_bystream(aap,domain,varargin{:},inpstream{tInd+1});
        
        % Other?
        otherimfn = '';
        if numel(inpstream) > (tInd+1)
            for s = 1:numel(inpstream)-(tInd+1)
                otherimfn = strvcat(otherimfn, aas_getfiles_bystream(aap,domain,varargin{:},inpstream{tInd+1+s}));
            end
        end
        
        %% Coregister
        % Coregister source to target
        x = spm_coreg(spm_vol(deblank(targetimfn)), spm_vol(sourceimfn), flags.estimate);
        % Set the new space for the mean EPI
        spm_get_space(sourceimfn, spm_matrix(x)\spm_get_space(sourceimfn));
                
        fprintf(['\t%s to %s realignment parameters:\n' ...
            '\tx: %0.4f   y: %0.4f   z: %0.4f   p: %0.4f   r: %0.4f   j: %0.4f\n'], ...
            inpstream{2}, inpstream{1}, ...
            x(1), x(2), x(3), x(4), x(5), x(6))

        %% Other
        if numel(inpstream) > (tInd+1)
            % Again, get space of source
            MM = spm_get_space(sourceimfn);
            
            for e = 1:size(otherimfn,1)
                % Apply the space of the coregistered source to the
                % remaining others (safest solution!)
                spm_get_space(deblank(otherimfn(e,:)), MM);
            end
        end
        
        %% Reslice:
        spm_reslice(strvcat(targetimfn,sourceimfn,otherimfn),flags.write);
        
        %% Describe the outputs and Diagnostics 
        % For source diag only
        if strcmp(aap.options.wheretoprocess,'localsingle')
           aas_checkreg(aap,domain,varargin{:},inpstream{tInd+1},inpstream{tInd});
        end
        
        % For outputs
        for s = 1:numel(outpstream)
            otherimfn = aas_getfiles_bystream(aap,domain,varargin{:},outpstream{s});
            otherimfn2 = '';
            for e = 1:size(otherimfn,1)
                [pth, fname, ext] = fileparts(deblank(otherimfn(e,:)));
                fpath = fullfile(pth, [flags.write.prefix fname ext]);

                % binarise if specified
                if isfield(aap.tasklist.currenttask.settings,'PVE') && ~isempty(aap.tasklist.currenttask.settings.PVE)
                    ninf = spm_vol(fpath);
                    Y = spm_read_vols(ninf);
                    Y = Y>=aap.tasklist.currenttask.settings.PVE;
                    nifti_write(fpath,Y,'Binarized',ninf)
                end
                
                otherimfn2 = strvcat(otherimfn2, fpath);
            end
            
            aap = aas_desc_outputs(aap,domain,varargin{:},outpstream{s},otherimfn2);
            if strcmp(aap.options.wheretoprocess,'localsingle')
                aas_checkreg(aap,domain,varargin{:},outpstream{s},inpstream{tInd});
            end
        end
        
    case 'checkrequirements'

end
