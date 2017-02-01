% FORMAT aa_export_toBIDS([aap,] dest[, 'anatt1', <stagetag[|streamname]>][, 'anatt2', <stagetag[|streamname]>])
% TODO: 
%   - multi-session
%   - multi-run
%   - covariates
%   - parametric
%   - multi-model

function aa_export_toBIDS(varargin)

if isstruct(varargin{1})
    aap = varargin{1};
    destpath = varargin{2};
    varargin(1:2) = [];
else
    loaded = load('aap_parameters');
    aap = loaded.aap;
    destpath = varargin{1};
    varargin(1) = [];
end
args = vargParser(varargin);

% source stages
stages = {aap.tasklist.main.module.name};
stage_anatt1 = cell_index(stages,'convert_structural') + cell_index(stages,'structuralfromnifti');
stage_anatt2 = cell_index(stages,'convert_t2');% + cell_index(stages,'t2fromnifti');
stage_func = cell_index(stages,'convert_epi') + cell_index(stages,'epifromnifti');
stage_fmap = cell_index(stages,'convert_fieldmap') + cell_index(stages,'fieldmapfromnifti');
stage_dwi = cell_index(stages,'convert_diffusion') + cell_index(stages,'diffusionfromnifti');

% folder
aas_makedir(aap,destpath);

% participant key file
fid_part = fopen(fullfile(destpath,'participants.tsv'),'w');
fprintf(fid_part,'participant_id'); % minimal
header_part = false; % header completed

% data
for subj = 1:numel(aap.acq_details.subjects)
    aas_log(aap,false,['INFO: Exporting subject: ' aas_getsubjname(aap,subj)])
    subjpath = fullfile(destpath,['sub-' aas_getsubjname(aap,subj)]);
    aas_makedir(aap,subjpath);
    if stage_anatt1
        aas_log(aap,false,'\tINFO: Exporting T1 anatomy')
        aap = aas_setcurrenttask(aap,stage_anatt1);
        aap.options.verbose = -1;
        fhdr = aas_getfiles_bystream(aap,'subject',subj,'structural_dicom_header','output');
        aap.options.verbose = 2;
        streamname = 'structural';
        if isfield(args,'anatt1'),
            out = textscan(args.anatt1,'%s','delimiter','|');
            currstage = out{1}{1};
            if numel(out{1}) > 1, streamname = out{1}{2}; end
            aap = aas_setcurrenttask(aap,aas_getmoduleindexfromtag(aap,currstage));
        end
        src = aas_getfiles_bystream(aap,'subject',subj,streamname,'output');
        
        % image
        if isempty(src)
            aas_log(aap,false,'WARNING: No image is available!')
            continue;
        end
        aas_makedir(aap,fullfile(subjpath,'anat')); 
        dest = fullfile(subjpath,'anat',sprintf('sub-%s_T1w.nii',aas_getsubjname(aap,subj)));
        copyfile(src,dest); gzip(dest); delete(dest);
        
        % header
        if isempty(fhdr)
            aas_log(aap,false,'WARNING: No header is available!')
        else
            loaded = load(fhdr); hdr = loaded.dcmhdr;
            json = struct(...
                'RepetitionTime',hdr.volumeTR,...
                'EchoTime',hdr.volumeTE,...
                'FlipAngle',hdr.FlipAngle...
                );
            savejson('',json,spm_file(dest,'ext','json'));
        end
    end
    if stage_anatt2
        aas_log(aap,false,'\tINFO: Exporting T2 anatomy')
        aap = aas_setcurrenttask(aap,stage_anatt2);
        aap.options.verbose = -1;
        fhdr = aas_getfiles_bystream(aap,'subject',subj,'t2_dicom_header','output');
        aap.options.verbose = 2;
        streamname = 't2';
        if isfield(args,'anatt2'),
            out = textscan(args.anatt2,'%s','delimiter','|');
            currstage = out{1}{1};
            if numel(out{1}) > 1, streamname = out{1}{2}; end
            aap = aas_setcurrenttask(aap,aas_getmoduleindexfromtag(aap,currstage));
        end        
        src = aas_getfiles_bystream(aap,'subject',subj,streamname,'output');
        
        % image
        if isempty(src)
            aas_log(aap,false,'WARNING: No image is available!')
            continue;
        end
        aas_makedir(aap,fullfile(subjpath,'anat')); 
        dest = fullfile(subjpath,'anat',sprintf('sub-%s_T2w.nii',aas_getsubjname(aap,subj)));
        copyfile(src,dest); gzip(dest); delete(dest);
        
        % header
        if isempty(fhdr)
            aas_log(aap,false,'WARNING: No header is available!')
        else
            loaded = load(fhdr); hdr = loaded.dcmhdr;
            if iscell(hdr), hdr = hdr{1}; end; hdr = hdr(1);
            json = struct(...
                'RepetitionTime',hdr.volumeTR,...
                'EchoTime',hdr.volumeTE,...
                'FlipAngle',hdr.FlipAngle...
                );
            savejson('',json,spm_file(dest,'ext','json'));
        end
    end
    if stage_func
        sliceaxes = {'ROW',{'i+' 'i-'}; 'COL',{'j+' 'j-'}};
        
        aap = aas_setcurrenttask(aap,stage_func);
        for sess = 1:numel(aap.acq_details.sessions)
            aas_log(aap,false,['\tINFO: Exporting fMRI session: ' aas_getsessname(aap,sess)])
            src = aas_getfiles_bystream(aap,'session',[subj sess],'epi','output');
            fhdr = aas_getfiles_bystream(aap,'session',[subj sess],'epi_dicom_header','output');
            
            % image
            if isempty(src)
                aas_log(aap,false,'WARNING: No image is available!')
                continue;
            end
            aas_makedir(aap,fullfile(subjpath,'func')); 
            dest = fullfile(subjpath,'func',sprintf('sub-%s_task-%s_bold.nii',aas_getsubjname(aap,subj),valueValidate(aap.acq_details.sessions(sess).name)));
            copyfile(src,dest); gzip(dest); delete(dest);
            
            % header
            loaded = load(fhdr); hdr = loaded.DICOMHEADERS{1};
            if isempty(fieldnames(hdr))
                aas_log(aap,false,'WARNING: No header information is available!')
            else
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
            end
            
            % events
            models = aap.tasksettings.aamod_firstlevel_model(1).model(2:end);
            selected_model = (strcmp({models.subject},aas_getsubjname(aap,subj)) | strcmp({models.subject},'*')) & ...
                (strcmp({models.session},aap.acq_details.sessions(sess).name) | strcmp({models.session},'*'));
            if ~any(selected_model)
                aas_log(aap,false,'WARNING: No model information is available!')
                continue
            end
            events = models(selected_model).event;
            ord = zeros(2,0);
            for e = 1:numel(events)
                ord = [ord vertcat(e*ones(1,numel(events(e).ons)),1:numel(events(e).ons))];
            end
            if size(events(1).ons,1) > size(events(1).ons,2) % column
                [ons,inde] = sort(vertcat(events.ons));
            else
                [ons,inde] = sort(horzcat(events.ons)');
            end
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
        
        aap = aas_setcurrenttask(aap,stage_fmap);
        for sess = 1:numel(aap.acq_details.sessions)
            aas_log(aap,false,['\tINFO: Exporting fieldmap for fMRI session: ' aas_getsessname(aap,sess)])
            fsrc = aas_getfiles_bystream(aap,'session',[subj sess],'fieldmap','output');
            aap.options.verbose = -1;
            fhdr = aas_getfiles_bystream(aap,'session',[subj sess],'fieldmap_dicom_header','output');
            aap.options.verbose = 2;
            
            % image
            if isempty(fsrc)
                aas_log(aap,false,'WARNING: No image is available!')
                continue;
            end
            aas_makedir(aap,fullfile(subjpath,'fmap')); 
            for f = 1:size(fsrc,1)
                src = deblank(fsrc(f,:));
                dest_format = fullfile(subjpath,'fmap',sprintf('sub-%s_run-%%d_%s.nii',aas_getsubjname(aap,subj),fmsuffix{f}));
                ind = 1;
                while exist(sprintf([dest_format '.gz'],ind),'file'), ind = ind + 1; end
                dest = sprintf(dest_format,ind);
                copyfile(src,dest); gzip(dest); delete(dest);
            end
            
            % header
            if isempty(fhdr)
                aas_log(aap,false,'WARNING: No header is available!')
            else
                loaded = load(fhdr); hdr = cell2mat(loaded.dcmhdr);
                TEs = unique([hdr.volumeTE]);
                if numel(TEs) < 2, TEs = sort(unique([hdr.EchoTime]),'ascend')/1000; end % try original (backward compatibility)
                json = struct(...
                    'RepetitionTime',hdr(1).volumeTR,...
                    'EchoTime1',TEs(1),...
                    'EchoTime2',TEs(2),...
                    'FlipAngle',hdr(1).FlipAngle,...
                    'IntendedFor',fullfile('func',sprintf('sub-%s_task-%s_bold.nii.gz',aas_getsubjname(aap,subj),valueValidate(aap.acq_details.sessions(sess).name))) ...
                    );
                savejson('',json,spm_file(dest,'ext','json'));
            end
        end
    end
    if stage_dwi
        sliceaxes = {'ROW',{'i+' 'i-'}; 'COL',{'j+' 'j-'}};
        
        aap = aas_setcurrenttask(aap,stage_dwi);
        for sess = 1:numel(aap.acq_details.diffusion_sessions)
            aas_log(aap,false,['\tINFO: Exporting DWI session: ' aas_getsessname(aap,sess)])
            src = aas_getfiles_bystream(aap,'diffusion_session',[subj sess],'diffusion_data','output');
            aap.options.verbose = -1;
            fhdr = aas_getfiles_bystream(aap,'diffusion_session',[subj sess],'diffusion_dicom_header','output');
            aap.options.verbose = 2;
            
            % image
            if isempty(src)
                aas_log(aap,false,'WARNING: No image is available!')
                continue;
            end
            aas_makedir(aap,fullfile(subjpath,'dwi'));
            dest = fullfile(subjpath,'dwi',sprintf('sub-%s_dwi.nii',aas_getsubjname(aap,subj)));
            copyfile(src,dest); gzip(dest); delete(dest);
            
            % header
            if isempty(fhdr)
                aas_log(aap,false,'WARNING: No header information is available!')
            else
                loaded = load(fhdr); hdr = loaded.DICOMHEADERS{1};
                json = struct(...
                    'RepetitionTime',hdr.volumeTR,...
                    'EchoTime',hdr.volumeTE,...
                    'FlipAngle',hdr.FlipAngle,...
                    'SliceTiming',hdr.slicetimes,...
                    'EffectiveEchoSpacing',hdr.echospacing,...
                    'PhaseEncodingDirection',[...
                    sliceaxes{cell_index(sliceaxes(:,1),deblank(hdr.InPlanePhaseEncodingDirection)),2}{aas_get_numaris4_numval(hdr.CSAImageHeaderInfo,'PhaseEncodingDirectionPositive')+1}...
                    ]...
                    );
                savejson('',json,spm_file(dest,'ext','json'));
            end
            % b
            copyfile(aas_getfiles_bystream(aap,'diffusion_session',[subj sess],'bvals','output'),...
                spm_file(dest,'ext','bval'));
            copyfile(aas_getfiles_bystream(aap,'diffusion_session',[subj sess],'bvecs','output'),...
                spm_file(dest,'ext','bvec'));
        end        
    end
    
    % participant key file
    hdr = hdr(1);
    entry = aas_getsubjname(aap,subj);
    
    if isfield(hdr,'PatientAge')
        if ~header_part, fprintf(fid_part,'\tage'); end
        out = textscan(hdr.PatientAge,'%03d%c'); [age, age_measure] = out{:};
        entry = sprintf('%s\t%d',entry,age);
    end
    if isfield(hdr,'PatientSex')
        if ~header_part, fprintf(fid_part,'\tsex'); end
        entry = sprintf('%s\t%c',entry,deblank(hdr.PatientSex));
    end
    
    if ~header_part, fprintf(fid_part,'\n'); header_part = true; end
    fprintf(fid_part,'%s\n',entry);
end

fclose(fid_part);

end

function out = valueValidate(in)
RULES = {...
    '-' '';
    '_' ''...
    };
out = strrep_multi(in,RULES(:,1),RULES(:,2));
end