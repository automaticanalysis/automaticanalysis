% FORMAT aa_export_toBIDS([aap,] dest)
% TODO: 
%   - diffusion
%   - multisession
%   - multirun
%   - covariates
%   - parametric

function aa_export_toBIDS(varargin)

if isstruct(varargin{1})
    aap = varargin{1};
    destpath = varargin{2};
else
    loaded = load('aap_parameters');
    aap = loaded.aap;
    destpath = varargin{1};
end

% source stages
stages = {aap.tasklist.main.module.name};
stage_anat = cell_index(stages,'convert_structural');
stage_func = cell_index(stages,'convert_epi');
stage_fmap = cell_index(stages,'convert_fieldmap');
stage_dwi = cell_index(stages,'convert_diffusion');

% data
aas_makedir(aap,destpath);
for subj = 1:numel(aap.acq_details.subjects)
    subjpath = fullfile(destpath,['sub-' aas_getsubjname(aap,subj)]);
    aas_makedir(aap,subjpath);
    if stage_anat
        aas_makedir(aap,fullfile(subjpath,'anat')); 
        aap = aas_setcurrenttask(aap,stage_anat);
        src = aas_getfiles_bystream(aap,'subject',subj,'structural','output');
        fhdr = aas_getfiles_bystream(aap,'subject',subj,'structural_dicom_header','output');
        
        % image
        dest = fullfile(subjpath,'anat',sprintf('sub-%s_T1w.nii',aas_getsubjname(aap,subj)));
        copyfile(src,dest); gzip(dest); delete(dest);
        
        % header
        loaded = load(fhdr); hdr = loaded.dcmhdr;
        json = struct(...
            'RepetitionTime',hdr.volumeTR,...
            'EchoTime',hdr.volumeTE,...
            'FlipAngle',hdr.FlipAngle...
            );
        savejson('',json,spm_file(dest,'ext','json'));
    end
    if stage_func
        sliceaxes = {'ROW',{'x+' 'x-'}; 'COL',{'y+' 'y-'}};
        
        aas_makedir(aap,fullfile(subjpath,'func')); 
        aap = aas_setcurrenttask(aap,stage_func);
        for sess = 1:numel(aap.acq_details.sessions)
            src = aas_getfiles_bystream(aap,'session',[subj sess],'epi','output');
            fhdr = aas_getfiles_bystream(aap,'session',[subj sess],'epi_dicom_header','output');
            
            % image
            dest = fullfile(subjpath,'func',sprintf('sub-%s_task-%s_bold.nii',aas_getsubjname(aap,subj),aap.acq_details.sessions(sess).name));
            copyfile(src,dest); gzip(dest); delete(dest);
            
            % header
            loaded = load(fhdr); hdr = loaded.DICOMHEADERS{1};
            json = struct(...
                'RepetitionTime',hdr.volumeTR,...
                'EchoTime',hdr.volumeTE,...
                'FlipAngle',hdr.FlipAngle,...
                'SliceTiming',hdr.slicetimes,...
                'EffectiveEchoSpacing',hdr.echospacing,...
                'PhaseEncodingDirection',[...
                    sliceaxes{cell_index(sliceaxes(:,1),deblank(hdr.InPlanePhaseEncodingDirection)),2}{aas_get_numaris4_numval(hdr.CSAImageHeaderInfo,'PhaseEncodingDirectionPositive')+1}...
                    ],...
                'TaskName',aap.acq_details.sessions(sess).name ...
                );
            savejson('',json,spm_file(dest,'ext','json'));
            
            % events
            models = aap.tasksettings.aamod_firstlevel_model.model(2:end);
            events = models((strcmp({models.subject},aas_getsubjname(aap,subj)) | strcmp({models.subject},'*')) & ...
                (strcmp({models.session},aap.acq_details.sessions(sess).name) | strcmp({models.session},'*'))).event;
            ord = zeros(2,0);
            for e = 1:numel(events)
                ord = [ord vertcat(e*ones(1,numel(events(e).ons)),1:numel(events(e).ons))];
            end
            [ons,inde] = sort([events.ons]);
            ord = ord(:,inde);
            lines = sprintf('onset\tduration\tweight\ttrial_type\n');
            for l = 1:size(ord,2)
                try dur = events(ord(1,l)).dur(ord(2,l)); catch, dur = events(ord(1,l)).dur; end
                lines = sprintf('%s%1.3f\t%1.3f\t%1.3f\t%s\n',lines,events(ord(1,l)).ons(ord(2,l)),dur,1,events(ord(1,l)).name);
            end
            fid = fopen(spm_file(strrep(dest,'_bold.','_events.'),'ext','tsv'),'w');
            fprintf(fid,lines);
            fclose(fid);
        end
    end
    if stage_fmap
        fmsuffix = {'magnitude1' 'magnitude2' 'phasediff' };
        
        aas_makedir(aap,fullfile(subjpath,'fmap')); 
        aap = aas_setcurrenttask(aap,stage_fmap);
        for sess = 1:numel(aap.acq_details.sessions)
            fsrc = aas_getfiles_bystream(aap,'session',[subj sess],'fieldmap','output');
            fhdr = aas_getfiles_bystream(aap,'session',[subj sess],'fieldmap_dicom_header','output');
            
            % image
            for f = 1:size(fsrc,1)
                src = deblank(fsrc(f,:));
                dest = fullfile(subjpath,'fmap',sprintf('sub-%s_task-%s_%s.nii',aas_getsubjname(aap,subj),aap.acq_details.sessions(sess).name,fmsuffix{f}));
                copyfile(src,dest); gzip(dest); delete(dest);
            end
            
            % header
            loaded = load(fhdr); hdr = cell2mat(loaded.dcmhdr);
            TEs = unique([hdr.volumeTE]);
            if numel(TEs) < 2, TEs = sort(unique([hdr.EchoTime]),'ascend')/1000; end % try original (backward compatibility)
            json = struct(...
                'RepetitionTime',hdr(1).volumeTR,...
                'EchoTime1',TEs(1),...
                'EchoTime2',TEs(2),...
                'FlipAngle',hdr(1).FlipAngle,...
                'IntendedFor',fullfile('func',sprintf('sub-%s_task-%s_bold.nii.gz',aas_getsubjname(aap,subj),aap.acq_details.sessions(sess).name)) ...
                );
            savejson('',json,spm_file(dest,'ext','json'));
        end
    end
end

end