% Automatic analysis - add subject to a planned analysis
% Usually in your user script.
% Using this function supercedes separately filling
% aap_acq_details.subjects and aap.acq_details_brukersessionnums. It is
% more convenient when you have many subjects as the correspondence between
% subject name and series numbers is more transparent.
%
%function [aap]=aas_addsubject(aap,name,seriesnumbers,ignoreseries,specialseries,args)
% name= subject filename (may include UNIX wildcards, e.g., CBU060500*/*)
% seriesnumbers=
%	DICOM source: series numbers of EPIs for this subject
%	NIFTI source: cell array containing (one or more)
%		- string: full or realtive path (from rawdatadir - only one is supported) to structural
%		- string: full or realtive path (from rawdatadir - only one is supported) to 4D NIFTI of one fMRI session
%		- string: full or realtive path (from rawdatadir - only one is supported) to wholebrain EPI (only after fMRI)
%		- cell array (nested): full path to 3D NIFTI files of one fMRI session
% ignoreseries parameter=series numbers of any series to be ignored in the
% analysis (e.g. a repeated structural) [added by djm 20/3/06]
% specialseries= special series to be converted
% args= 1x2N cell, containing N arguments (odd: name, even: val). E.g.:
% 	{'--structural', serienum}: specifying serienumber for structural

function [aap]=aas_addsubject(aap,name,seriesnumbers,ignoreseries,specialseries,diffusion_seriesnumbers,args)

% Blank template for a subject entry
f=fieldnames(aap.schema.acq_details.subjects);
for field=1:length(f)
    if (~strcmp(f{field},'ATTRIBUTE'))
        thissubj.(f{field})=[];
    end;
end;

% Both MEG & MRI or just MRI?
try
    if (length(name)==2)
        thissubj.megname=name{1};
        thissubj.mriname=name{2};
    else
        thissubj.mriname=name;
    end;
catch
    aas_log(aap,true,'In aas_addsubject, expecting either single name for MRI in single quotes, or two names for MEG written like this {''megname'',''mriname''}.');
end;
try
    if isnumeric(seriesnumbers) || isnumeric(seriesnumbers{1}) % DICOM series number
        thissubj.seriesnumbers=seriesnumbers;
    else
        fMRI = {}; MEG = {};
        fMRIdim = [];
        for s = 1:numel(seriesnumbers)
            if iscell(seriesnumbers{s}) % multiple 3D files
                fMRI{end+1} = seriesnumbers{s};
            elseif ischar(seriesnumbers{s}) % single NIFTI file
                if exist(seriesnumbers{s},'file') % full path
                    V = spm_vol(seriesnumbers{s});
                else
                    if ~isempty(thissubj.megname) % try meg
                        if exist(fullfile(meg_findvol(aap,thissubj.megname,1),seriesnumbers{s}),'file')
                            MEG{end+1} = seriesnumbers{s};
                            continue;
                        end
                    end
                    % try to find in rawdatadir
                    V = spm_vol(fullfile(aas_findvol(aap,''),seriesnumbers{s}));
                end
                if numel(V) > 1 % 4D --> fMRI
                    fMRI{end+1} = seriesnumbers{s};
                    fMRIdim = V(1).dim;
                else % 3D --> structural
                    if ~isempty(fMRIdim) && all(fMRIdim(1:2) == V.dim(1:2)) % same inplane resolution as fMRI
                        thissubj.wholebrain_epi=seriesnumbers(s);
                    else
                        thissubj.structural=seriesnumbers(s);
                    end
                end
            else isnumeric(seriesnumbers{s}) % mixed: DICOM series number for fMRI
                thissubj.seriesnumbers=seriesnumbers{s};
            end
        end
        if ~isempty(fMRI)
            thissubj.seriesnumbers=fMRI;
        end
        if ~isempty(MEG)
            thissubj.megseriesnumbers=MEG;
        end
    end
catch
end;


% [djm 20/3/06]
if nargin>=4
    thissubj.ignoreseries=ignoreseries;
end;

if nargin>=5
    thissubj.specialseries=specialseries;
end;

if nargin>=6
    thissubj.diffusion_seriesnumbers=diffusion_seriesnumbers;
end;

if nargin >=7 % extra parameters
    structural = get_arg(args,'--structural');
    if ~isempty(structural)
        aap.options.autoidentifystructural_choose.subject(end+1) = struct(...
            'name',mri_findvol(aap,name),...
            'serie',structural);
    end
end

% And put into acq_details, replacing a single blank entry if it exists
if (length(aap.acq_details.subjects)==1 && isempty(aap.acq_details.subjects.mriname) && isempty(aap.acq_details.subjects.megname))
    aap.acq_details.subjects=thissubj;
else
    aap.acq_details.subjects(end+1)=thissubj;
end;
end

function val = get_arg(args,name)
    val = [];
    ind = find(strcmp(args,name));
    if ~isempty(ind)
        val = args{ind+1};
    end
end

