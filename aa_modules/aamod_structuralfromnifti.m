% AA module - structural from NIFTI

function [aap,resp]=aamod_structuralfromnifti(aap,task,subj)

resp='';

switch task
    case 'report'
    case 'doit'
        stream = aap.tasklist.currenttask.outputstreams.stream{1};
        
        % Select
        series = horzcat(aap.acq_details.subjects(subj).(stream){:});
        if ~iscell(series) || ~ischar(series{1})
            aas_log(aap,true,'ERROR: Was expecting list of filename(s) in cell array');
        end
        switch numel(series)
            case 1
                series = series{1};
            otherwise
                if aap.options.autoidentifystructural_choosefirst, series = series{1}; end
                if aap.options.autoidentifystructural_chooselast, series = series{end}; end
        end
        
        % Process
        niftifile = series;
        if ~exist(niftifile,'file')
            niftisearchpth=aas_findvol(aap,'');
            if ~isempty(niftisearchpth)
                niftifile = fullfile(niftisearchpth,series); % only the first file is processed
            end
        end
        comp = false;
        if strcmp(spm_file(niftifile,'Ext'),'gz'),
            comp = true;
            gunzip(niftifile);
            niftifile = niftifile(1:end-3);
        end
        V=spm_vol(niftifile);
        sesspth=fullfile(aas_getsubjpath(aap,subj),aap.directory_conventions.structdirname);
        aas_makedir(aap,sesspth);
        [pth fle ext]=fileparts(series);
        Y=spm_read_vols(V(1));
        fn=fullfile(sesspth,[sprintf('%s_%04d',fle,1) '.nii']);
        V(1).fname=fn;
        V(1).n=[1 1];
        % Write out the files, now likely in 3d
        if comp, delete(niftifile); end
        spm_write_vol(V(1),Y);        
        aap=aas_desc_outputs(aap,subj,stream,fn);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end