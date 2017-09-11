function aap = aas_addsubject(aap, varargin)
% Add subject/data to the analysis. It can be called it multiple times to add more subjects and/or more sources to a particular subject.
%
% FORMAT function aap = aas_addsubject(aap, data)
% Process only autoidentified images. User will be a warned about the lack of series specification. Subject name, which is used as a reference in aa 
% (i.e. when calling aas_addevent, aas_addcontrasts, etc.), will be automatically defined based on aap.directory_conventions.subject_directory_format:
%   0   - based on predefined list stored in aap.directory_conventions.subject_directory_names
%   1   - based on the "data" (see below)
%   2   - based on the order of specification (S01, S02, etc.)
%
% aap           - aap structure with parameters and tasklist
% data          - subject foldername within database. It should correspond to aap.directory_conventions.subjectoutputformat or 
%                 aap.directory_conventions.megsubjectoutputformat
%
%
% FORMAT function aap = aas_addsubject(aap, name, data)
% As above, but manually specifies subject name. aap.directory_conventions.subject_directory_format must be set to 3.
%
% name          - subject name
%
%
% FORMAT function aap = aas_addsubject(___,'name',subjectname)
% Another way to specify subject name, which ignores aap.directory_conventions.subject_directory_format.
%
% subjectname   - subject name
%
%
% FORMAT function aap = aas_addsubject(___,'functional',series)
% Specify functional (fMRI/MEG) data.
%
% series        - for DICOM: array of series number(s) of EPIs. E.g.: 
%                   two series of single-echo EPI: [5 10]
%                   two series of single-echo EPI and one series of multi-echo EPI with 5 echos: {5 10 15:19}
%               - for NIfTI: cell array containing one or more
%                   for structural: string of full or realtive path (from a single rawdatadir)
%                   for 4D NIfTI: string of full or realtive path (from a single rawdatadir)
%                   for wholebrain EPI: string of full or realtive path (from a single rawdatadir). Can be specified only after fMRI series.
%                   for 3D NIfTI: cell array (i.e. nested) of strings of full path
%                 Strings can be replaced by structures with fields 'fname' (path to image) and 'hdr' (path to header) to specify metadata.
%               - for MEG: filename of the *.fif file in the subject folder (i.e. no full path).
% Missing series can be specified either with "0" (for numerical array input) or with "[]" (for cell array input).
%
%
% FORMAT function aap = aas_addsubject(___,'diffusion',series)
% Specify diffusion-weighted MRI data.
%
% series        - for DICOM: numeric array of series number(s)
%               - for NIfTI: cell of structure(s) with fields 'fname' (path to image), and 'bval', 'bvec'(path to bvals and bvecs)
%
%
% FORMAT function aap = aas_addsubject(___,'structural', series)
% Specify structural data (overwrites autoidentification).
%
% series        - for DICOM: numeric array of series number
%               - for NIfTI: cell of string (path to image) or cell of structure with fields 'fname' (path to image) and 'hdr' (path to header)
%
%
% FORMAT function aap = aas_addsubject(___,'fieldmaps', series)
% Specify fieldmap data (overwrites autoidentification).
%
% series        - for DICOM: numeric array of series numbers
%               - for NIfTI: cell of structure with fields 'fname' (cell of 3 filenames - 2x magnitude + 1x phase) and 'hdr' (path to header)
%
%
% FORMAT function aap = aas_addsubject(___,'specialseries', series)
% Specify 'special' data (e.g. ASL, MTI, MPM).
%
% series        - for DICOM: cell array of numeric arrays of series numbers
%               - for NIfTI: not supported yet
%
%
% FORMAT function aap = aas_addsubject(___,'ignoreseries', series)
% Specify DICOM series to be ignored during autoidentification.
%
% series        - numeric arrays of series numbers

%% Parse
iMRIData = 1; % new subject
iMEGData = 1;
if isempty(aap.acq_details.subjects(end).subjname)
    subjind = 1;
else
    subjind = numel(aap.acq_details.subjects) + 1;
end

name = '';
switch aap.directory_conventions.subject_directory_format
    case 0 % from predefined list
        name = aap.directory_conventions.subject_directory_names{subjind};
        data = varargin{1};
        varargin(1)= [];
    case 1 % from data
        data = varargin{1};
        varargin(1)= [];
    case 2 % S#
        name = sprintf('S%02d',subjind);
        data = varargin{1};
        varargin(1)= [];
    case 3 % manual
        name = varargin{1};
        data = varargin{2};
        varargin(1:2)= [];
    otherwise
        aas_log(aap,true,sprintf('ERROR: Unknown subject directory format (aap.directory_conventions.subject_directory_format=%d. Value only 0-3 is allowed.',aap.directory_conventions.subject_directory_format));
end
try
    args = vargParser(varargin);
catch
    aas_log(aap,false,sprintf('ERROR in %s: incorrect arguments',mfilename),'Errors')
    help(mfilename);
    error('ERROR in %s: incorrect arguments',mfilename);
end

%% Sanity, compatiblity check
if isempty(varargin)
    aas_log(aap,false,'WARNING: No series has been specified!\n')
else
    if ~isa(varargin{1},'char')
        aas_log(aap,true,sprintf('ERROR: Arguments are  different from what expected!\n %s',help('aas_addsubject')))
    end
end

%% Initialize subject
% with a blank template for a subject entry
fields=fieldnames(aap.schema.acq_details.subjects);
fields(strcmp(fields,'ATTRIBUTE')) = [];
for field=fields'
    thissubj.(field{1})={[]};
end
fields(strcmp(fields,'subjname')) = [];

% search for existing subject
if isfield(args,'name'), name = args.name; end
if ~isempty(name) && ~isempty(aap.acq_details.subjects(1).subjname)
% name specified --> check whether subject already exists (only if there is at least one already)    
    subjserach = cell_index({aap.acq_details.subjects.subjname},name);
    if subjserach
        subjind = subjserach; 
        thissubj = aap.acq_details.subjects(subjind);
        iMRIData = numel(thissubj.mriname)+1;
        iMEGData = numel(thissubj.megname)+1;
        for field=fields'
            thissubj.(field{1}){end+1}=[];
        end
    end
end

%% Data
try
    if iscell(data) && numel(data) == 2 % MEG
        thissubj.megname{iMEGData}=data{1};
        thissubj.mriname{iMRIData}=data{2};
        if isempty(name), name = aas_megname2subjname(aap,sprintf(aap.directory_conventions.megsubjectoutputformat,thissubj.megname{1})); end
    else % MRI
        thissubj.mriname{iMRIData}=data;
        if isempty(name), name = aas_mriname2subjname(aap,sprintf(aap.directory_conventions.subjectoutputformat,thissubj.mriname{1})); end
    end
catch
    aas_log(aap,true,'In aas_addsubject, name is expected to be either single string for MRI, or a cell of two for MEG written like this {''megname'',''mriname''}.');
end

thissubj.subjname = name;

%% Series
if isfield(args,'functional') && ~isempty(args.functional) 
    if isnumeric(args.functional) || isnumeric(args.functional{1}) % DICOM series number --> MRI
        thissubj.seriesnumbers{iMRIData}=args.functional;
    else
        fMRI = {}; MEG = {};
        fMRIdim = [];
        for s = 1:numel(args.functional)
            if iscell(args.functional{s}) % multiple 3D files
                fMRI{end+1} = args.functional{s};
            elseif ischar(args.functional{s}) ||... % NIfTI file
                    isstruct(args.functional{s})    % hdr+fname
                % Get filename
                
                if isstruct(args.functional{s})
                    if numel(args.functional{s}.fname) > 1 % multiple 3D files
                        fMRI{end+1} = args.functional{s};
                        continue;
                    end
                    fname = args.functional{s}.fname;
                else
                    fname = args.functional{s};
                end
                
                % - try in rawmegdatadir
                if ~exist(fname,'file')
                    if ~isempty(thissubj.megname{iMEGData})
                        tmpaap = aap;
                        tmpaap.directory_conventions.megsubjectoutputformat = '%s';
                        if exist(fullfile(meg_findvol(aap,thissubj.megname{iMEGData},'fp'),fname),'file') ||...
                                ~isempty(meg_findvol(tmpaap,fname)) % try empty room
                            MEG{end+1} = fname;
                            continue;
                        end
                    end
                end
                % - try in rawdatadir
                if ~exist(fname,'file'), fname = fullfile(aas_findvol(aap,spm_file(fname,'path')),spm_file(fname,'filename')); end
                
                if ~exist(fname,'file'), aas_log(aap,1,sprintf('ERROR: File %s does not exist!',fname)); end
                
                % Sort
                if strcmp(spm_file(fname,'ext'),'fif')
                     MEG{end+1} = fname;
                     continue;
                end
                
                V = spm_vol(fname);
                if numel(V) > 1 % 4D --> fMRI
                    fMRI{end+1} = args.functional{s};
                    fMRIdim = V(1).dim;
                else % 3D --> structural
                    if ~isempty(fMRIdim) && all(fMRIdim(1:2) == V.dim(1:2)) % same inplane resolution as fMRI
                        thissubj.wholebrain_epi{iMRIData}=args.functional(s);
                    else
                        thissubj.structural{iMRIData}=args.functional(s);
                    end
                end
            elseif isempty(args.functional{s}) % missing series
                fMRI{end+1} = [];
            else % mixed: DICOM series number for fMRI
                thissubj.seriesnumbers{iMRIData}=args.functional{s};
            end
        end
        if ~isempty(fMRI)
            thissubj.seriesnumbers{iMRIData}=fMRI;
        end
        if ~isempty(MEG)
            thissubj.megseriesnumbers{iMEGData}=MEG;
        end
    end
end

if isfield(args,'diffusion') && ~isempty(args.diffusion)
    thissubj.diffusion_seriesnumbers{iMRIData}=args.diffusion;
end

for meas = {'structural' 'fieldmaps' 'specialseries' 'ignoreseries'}
    if isfield(args,meas{1}) && ~isempty(args.(meas{1}))
        thissubj.(meas{1}){iMRIData}=args.(meas{1});
    end
end

% And put into acq_details, replacing a single blank entry if it exists
aap.acq_details.subjects(subjind)=thissubj;
end