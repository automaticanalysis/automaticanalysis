% Automatic analysis - add subject to the analysis. 
% You may call it multiple times to add more (DICOM) sources to a particular subject (MRI only!).
%
% FORMAT function aap = aas_addsubject(aap, [name], data, [['name',<subject name>], ['functional',<functional series>], ['diffusion',<diffusion series>], ['structural',<structural series>], ['fieldmaps',<fieldmaps series>], ['specialseries',<specialseries>], ['ignoreseries',<ignoreseries>]])
%   - aap: aap structure with parameters and tasklist
%   - name: subject name used as a reference in aa (i.e. when calling aas_addevent, aas_addcontrasts, etc.) | aap.directory_conventions.subject_directory_format == 3
%       When not specified, it will be automatically defined 
%         - based on predefined list stored in aap.directory_conventions.subject_directory_names | aap.directory_conventions.subject_directory_format == 0
%         - based on the "data" | aap.directory_conventions.subject_directory_format == 1
%         - based on subject number (S01, S02, etc.) | aap.directory_conventions.subject_directory_format == 2
%   - data: subject foldername within database (may include UNIX wildcards, e.g., CBU060500*/*)
%
%   'name', subjname: manual entry of subject name (overwrites aap.directory_conventions.subject_directory_format)
%   'functional', seriesnum: functional MRI or MEG series
%     For fMRI, seriesnum can be:
%         DICOM source: series numbers of EPIs for this subject
%         NIFTI source: cell array containing (one or more)
%           - string: full or realtive path (from rawdatadir - only one is supported) to structural
%           - string: full or realtive path (from rawdatadir - only one is supported) to 4D NIFTI of one fMRI session
%           - string: full or realtive path (from rawdatadir - only one is supported) to wholebrain EPI (only after fMRI)
%           - cell array (nested): full path to 3D NIFTI files of one fMRI session
%         all strings can be structures with fields 'fname' (path to image), and 'hdr', 'bval', 'bvec' (path to header, bvals and bvecs)
%     For MEG, seriesnum is the filename of the *.fif file.
%   'diffusion', seriesnum: diffusion-weighted MRI series
%   'structural', serienum or cellstr of 3D NIFTI filename: specifying structural
%   'fieldmaps', serienums or cell of structure with fields 'fname' (3 (2x mag + 1x phase)x 3D NIFTI filenames) and 'hdr' (path to header): specifying fieldmap
%   'specialseries', seriesnum: special series to be converted
%   'ignoreseries', seriesnum: series to be ignored in the analysis (e.g. a repeated structural) [added by djm 20/3/06]

function aap = aas_addsubject(aap, varargin)

%% Parse
iMRIData = 1; % new subject
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
    return
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
fields(strcmp(fields,'megname')) = [];
thissubj.megname = '';

% search for existing subject
if isfield(args,'name'), name = args.name; end
if ~isempty(name) % name specified --> check whether subject already exists
    subjserach = cell_index({aap.acq_details.subjects.subjname},name);
    if subjserach
        subjind = subjserach; 
        thissubj = aap.acq_details.subjects(subjind);
        iMRIData = numel(thissubj.mriname)+1;
        for field=fields'
            thissubj.(field{1}){end+1}=[];
        end
    end
end

%% Data
try
    if iscell(data) && numel(data) == 2 % MEG
        thissubj.megname=data{1};
        thissubj.mriname{iMRIData}=data{2};
        if isempty(name), name = aas_megname2subjname(aap,sprintf(aap.directory_conventions.megsubjectoutputformat,thissubj.megname)); end
    else % MRI
        thissubj.mriname{iMRIData}=data;
        if isempty(name), name = aas_mriname2subjname(aap,sprintf(aap.directory_conventions.subjectoutputformat,thissubj.mriname{1})); end
    end
catch
    aas_log(aap,true,'In aas_addsubject, name is expected to be either single string for MRI, or a cell of two for MEG written like this {''megname'',''mriname''}.');
end

thissubj.subjname = name;

%% Series
if isfield(args,'functional')
    if isnumeric(args.functional) || isnumeric(args.functional{1}) % DICOM series number
        thissubj.seriesnumbers{iMRIData}=args.functional;
    else
        fMRI = {}; MEG = {};
        fMRIdim = [];
        for s = 1:numel(args.functional)
            if iscell(args.functional{s}) % multiple 3D files
                fMRI{end+1} = args.functional{s};
            elseif ischar(args.functional{s}) || isstruct(args.functional{s}) % single NIFTI file or hdr+fname
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
                
                % Check if exists full path
                % - try meg
                if ~exist(fname,'file')
                    if ~isempty(thissubj.megname)
                        tmpaap = aap;
                        tmpaap.directory_conventions.megsubjectoutputformat = '%s';
                        if exist(fullfile(meg_findvol(aap,thissubj.megname,'fp'),fname),'file') ||...
                                ~isempty(meg_findvol(tmpaap,fname)) % try empty room
                            MEG{end+1} = fname;
                            continue;
                        end
                    end
                end
                % - try to find in rawdatadir
                if ~exist(fname,'file'), fname = fullfile(aas_findvol(aap,''),fname); end
                if ~exist(fname,'file'), aas_log(aap,1,sprintf('ERROR: File %s does not exist!',fname)); end
                
                % Sort
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
            else isnumeric(args.functional{s}) % mixed: DICOM series number for fMRI
                thissubj.seriesnumbers{iMRIData}=args.functional{s};
            end
        end
        if ~isempty(fMRI)
            thissubj.seriesnumbers{iMRIData}=fMRI;
        end
        if ~isempty(MEG)
            thissubj.megseriesnumbers=MEG;
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