function aa_export_toBIDS(varargin)
%
% BIDS export -- create a BIDS-compliant directory tree from an extant aa analysis
%
% This function can run standalone or be included at the end of a userscript.
%
% usage: aa_export_toBIDS([aap,] dest[, 'anatt1', <stagetag[|streamname]>][, 'anatt2', <stagetag[|streamname]>][, 'model', <stageindex>|'session'])
%
% also uses: aap.directory_conventions.BIDSfiles
%
% example usage:
%
%	1) aa_export_toBIDS('/put/results/here');
%	2) aa_export_toBIDS(aap,'/put/results/here');
%	3) aa_export_toBIDS('put/results/here','anatt1','aamod_freesurfer_deface_00001|defaced_structural')
%	4) aa_export_toBIDS('put/results/here','model',1)
%	5) aa_export_toBIDS('put/results/here','model','session')
%	
% notes
%
% a) if you pass the aap struct, you must first load it from aap_parameters.mat
% (this version contains additional required fields). If you don't pass aap, it 
% will be loaded from aap_parameters.m for you. As such, there really is no 
% reason you would ever pass in aap, assuming you first cd to the directory
% where aap_parameters.mat lives
%
% b) example usage #3 shows how to override the default structural
% (i.e. t1) file to copy (in this example a defaced structural is used, which
% is probably required if you're going to upload results to a public
% repository like openfMRI. The syntax goes: stage|streamname. If you don't
% supply an explict stream for t1 (or, optionally, t2), the first appearance
% of t1 [t2] in the analysis tasklist will be used, which is probably the
% output from dicom covert (which is not defaced).
%
% c) eample usage #4 and #5 show how to select model to extract events. aa
% uses the first model by default. A number here can specify a certain
% model, while 'session' tells aa that every session has its own specific model.
%
% d) BIDS requires three files to appear in the top level directory:
%
%		README - a plaintext (ASCII or UTF-8) description of the data
%		CHANGES- a plaintext (ASCII or UTF-8) list of version changes
%		dataset_description.json - a JSON description of the data (see the
%			current specification for required and optional fields)
%
%	You must provide these files to be BIDS-compliant. This function will
%	attempt to copy them from aap.directory_conventions.BIDSfiles, if the
%	field is defined and the directory exists (otherwise you'll have to add 
%	them by hand). Note there are a number of optional files that can be also
%	be included at the top level -- for convenience, all files that live in 
%	aap.directory_convention.BIDSfiles will be copied for you.
%
% See bids.neuroimaging.io/bids_spec1.0.0.pdf for specification details.
%
% This function is a work in progress. TODO:
%
%   - multi-model
%   - multi-session
%   - covariates
%   - parametric

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
if ~isfield(args,'model'), args.model = 1; end % default model

% identify source stages present from the tasklist

stages = {aap.tasklist.main.module.name};

stage_anatt1 = cell_index(stages,'convert_structural') + cell_index(stages,'structuralfromnifti');
stage_anatt2 = cell_index(stages,'convert_t2');% + cell_index(stages,'t2fromnifti');
stage_func = cell_index(stages,'convert_epi') + cell_index(stages,'epifromnifti');
stage_fmap = cell_index(stages,'convert_fieldmap') + cell_index(stages,'fieldmapfromnifti');
stage_dwi = cell_index(stages,'convert_diffusion') + cell_index(stages,'diffusionfromnifti');

% create the top-level directory

aas_makedir(aap, destpath);

% copy required BIDS files if they exist

success = 0;
if isfield(aap.directory_conventions,'BIDSfiles') && ~isempty(aap.directory_conventions.BIDSfiles)
	bidsdir = fullfile(aap.acq_details.root, aap.directory_conventions.BIDSfiles);
	if exist(bidsdir,'dir')
		success = ~system(sprintf('cp %s/* %s', bidsdir, destpath));
		% these three files are compulsory
		if ~exist(fullfile(destpath,'dataset_description.json'),'file'); success=0; end
		if ~exist(fullfile(destpath,'README'),'file'); success=0; end
		if ~exist(fullfile(destpath,'CHANGES'),'file'); success=0; end
	end
end

if ~success
	aas_log(aap, false, 'WARNING: Required top-level files dataset_description.json, README, and/or CHANGES not copied -- you must add these by hand.');
end

% participant key file

fid_part = fopen(fullfile(destpath,'participants.tsv'),'w');
fprintf(fid_part,'participant_id'); % minimal
header_part = false; % header completed

% data

for subj = 1:numel(aap.acq_details.subjects)
	
    aas_log(aap,false,['INFO: Exporting subject: ' aas_getsubjname(aap,subj)])
    subname = aas_getsubjname(aap,subj);
    suboutname = subname;
    if ~isequal(suboutname(1:4),'sub-')
        suboutname = ['sub-' suboutname];
	end
	
    subjpath = fullfile(destpath,suboutname);
	
    aas_makedir(aap,subjpath);
	
    if stage_anatt1
		
        aas_log(aap,false,'\tINFO: Exporting T1 anatomy')
        aap = aas_setcurrenttask(aap,stage_anatt1);
        aap.options.verbose = -1;
        fhdr = aas_getfiles_bystream(aap,'subject',subj,'structural_dicom_header','output');
        aap.options.verbose = 2;
        streamname = 'structural';

		if isfield(args,'anatt1')
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
        dest = fullfile(subjpath,'anat',sprintf('%s_T1w.nii',suboutname));
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
	
    if stage_anatt2
		
        aas_log(aap,false,'\tINFO: Exporting T2 anatomy')
        aap = aas_setcurrenttask(aap,stage_anatt2);
        aap.options.verbose = -1;
        fhdr = aas_getfiles_bystream(aap,'subject',subj,'t2_dicom_header','output');
        aap.options.verbose = 2;
        streamname = 't2';
        if isfield(args,'anatt2')
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
        dest = fullfile(subjpath,'anat',sprintf('%s_T2w.nii',suboutname));
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
            taskname = aas_getsessname(aap,sess);
            aas_log(aap,false,['\tINFO: Exporting fMRI session: ' taskname])
            src_main = aas_getfiles_bystream(aap,'session',[subj sess],'epi','output');

            % image
            if isempty(src_main)
                aas_log(aap,false,'WARNING: No image is available!')
                continue;
            end
            
            src = src_main;
            if any(strcmp(aas_getstreams(aap,'output'),'dummyscans'))
                src_dummy = aas_getfiles_bystream(aap,'session',[subj sess],'dummyscans','output');
                if ~isempty(src_dummy)
                    src = spm_file(tempname,'ext','nii');
                    spm_file_merge(char(src_dummy,src_main),src);
                end
            end
            fhdr = aas_getfiles_bystream(aap,'session',[subj sess],'epi_dicom_header','output');
            
            isRun = regexp(taskname,'[_-]?[rR]un[0-9]*$');
            if isRun
                taskname = taskname(1:isRun-1);
                runNo = str2double(regexp(aap.acq_details.sessions(sess).name(isRun:end),'[0-9]*','match'));
            end
            
            aas_makedir(aap,fullfile(subjpath,'func')); 
            dest = fullfile(subjpath,'func',sprintf('%s_task-%s_bold.nii',suboutname,valueValidate(taskname)));
            if isRun, dest = strrep(dest,'bold.nii',sprintf('run-%d_bold.nii',runNo)); end
            copyfile(src,dest); gzip(dest); delete(dest); if (exist('src_dummy','var') && ~isempty(src_dummy)), delete(src); clear src_dummy; end       
            % header
            loaded = load(fhdr); hdr = loaded.DICOMHEADERS{1};
            if isempty(fieldnames(hdr))
                aas_log(aap,false,'WARNING: No header information is available!')
            else
            if iscell(hdr), hdr = hdr{1}; end; hdr = hdr(1);
                json = struct(...
                    'RepetitionTime',hdr.volumeTR,...
                    'EchoTime',hdr.volumeTE,...
                    'FlipAngle',hdr.FlipAngle,...
                    'SliceTiming',hdr.slicetimes,...
                    'EffectiveEchoSpacing',hdr.echospacing,...
                    'PhaseEncodingDirection',[...
                    sliceaxes{cell_index(sliceaxes(:,1),deblank(hdr.InPlanePhaseEncodingDirection)),2}{aas_get_numaris4_numval(hdr.CSAImageHeaderInfo,'PhaseEncodingDirectionPositive')+1}...
                    ],...
                    'TaskName',taskname ...
                    );
                savejson('',json,spm_file(dest,'ext','json'));
            end
            
            % events
			
            stageModels = aap.tasklist.main.module(strcmp({aap.tasklist.main.module.name},'aamod_firstlevel_model'));
            
			if isempty(stageModels)
				aas_log(aap,false,'WARNING: No model information is available!')
				continue
            end

            switch args.model
                case 'session'
                    stageindex = strcmp(arrayfun(@(x) x.extraparameters.aap.acq_details.selected_sessions, stageModels,'UniformOutput',false),aas_getsessname(aap,sess));
                otherwise
                    stageindex = args.model;
            end
            
            models = aap.tasksettings.aamod_firstlevel_model(stageindex).model(2:end);
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

            % trying to determine if events are a row or col vect
            % this can fail if events.ons contains one element...
            
            try
                [ons,inde] = sort(vertcat(events.ons));
            catch
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
                dest_format = fullfile(subjpath,'fmap',sprintf('%s_run-%%02d_%s.nii',suboutname,fmsuffix{f}));
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
                    'IntendedFor',fullfile('func',sprintf('%s_task-%s_bold.nii.gz',suboutname,valueValidate(aap.acq_details.sessions(sess).name))) ...
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
            src_main = aas_getfiles_bystream(aap,'diffusion_session',[subj sess],'diffusion_data','output');

            % image
            if isempty(src_main)
                aas_log(aap,false,'WARNING: No image is available!')
                continue;
            end
            
            if any(strcmp(aas_getstreams(aap,'output'),'dummyscans'))
                src = spm_file(tempname,'ext','nii');
                src_dummy = aas_getfiles_bystream(aap,'session',[subj sess],'dummyscans','output');
                spm_file_merge(char(src_dummy,src_main),src);
            else
                src = src_main;
            end
            aap.options.verbose = -1;
            fhdr = aas_getfiles_bystream(aap,'diffusion_session',[subj sess],'diffusion_dicom_header','output');
            aap.options.verbose = 2;
            
            aas_makedir(aap,fullfile(subjpath,'dwi'));
            dest = fullfile(subjpath,'dwi',sprintf('%s_dwi.nii',suboutname));
            copyfile(src,dest); gzip(dest); delete(dest); if exist('src_dummy','var'), delete(src); clear src_dummy; end
            
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
    entry = suboutname;
    
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

