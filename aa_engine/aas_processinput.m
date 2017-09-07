% Allow specify Subjects, Sessions and Events based on text file.
% Required options in aap.acq_details.input (see also aap_parameters_defaults.xml):
%	- list: text file:
%		Required format:
%		- a CSV file: cells are separated with semicolon
%		- subcells are separated with "_"
%		- a header in the first line
%	 	Required columns (may contain more): 
%		- "ID": subject identifiers
%		- "FMRI1": cells for the fMRI measurement; it must contain the same number of subcells as the coresponding column header
%			- first subcell is the CBU volunteer number (without "CBU")
%			- series numbers are defined in the consecutive subcells with the same order as their names in the coresponding column header
%		E.g.:
%			ID;Age;Sex;FMRI1_Loc_Con_LDeasy_LDhard
%			02;29;f;90952_7_6_4_5
%
%	- selected_sessions: subselection of a subset of series:
%		E.g. (based on the example above): [3 4] --> LDeasy and LDhard
%
%	- referencedirectory_tmpl: a template for path to the folder containing Event files (.mat files in SPM-format)
%		Required format:
%		- must contain "*" which will be substituted with "ID" (see above)
%		E.g.: '/imaging/ta02/ActionWords/Analysis/E-Prime/DAW*/ref'
%
% Tibor Auer MRC CBU Cambridge 2012-2013

function aap = aas_processinput(aap)

if ~aap.directory_conventions.continueanalysis
    sfx = 1;
    analysisid0 = aap.directory_conventions.analysisid;
    while exist(fullfile(aap.acq_details.root,aap.directory_conventions.analysisid),'dir')
        sfx = sfx + 1;
        aap.directory_conventions.analysisid = [analysisid0 '_' num2str(sfx)];
    end
end

% get module name for first-level analysis
modulenames = fieldnames(aap.tasksettings);
stagenumModel = cell_index(modulenames,'firstlevel_model');
if any(stagenumModel), firstlevel = [modulenames{stagenumModel} '_*']; end

% Correct EV onstets for number of dummies?
numdummies = aap.acq_details.input.correctEVfordummies*aap.acq_details.numdummies; 

% Reads in header
head = study_header(aap.acq_details.input.list);
LIST = importdata(aap.acq_details.input.list);
for v = 2:size(LIST,1)
    ID = list_index(LIST{v},head.ID);
    VOL = list_index(LIST{v},head.FMRI1);
    if numel(regexp(VOL,'[0-9]')) == numel(VOL), VOL = str2double(VOL); end
    nSess = 2; aSess = [];
    while ~isempty(list_index(LIST{1},head.FMRI1,nSess))
        aSess = horzcat(aSess,str2double(list_index(LIST{v},head.FMRI1,nSess)));
        nSess = nSess + 1;
    end
    nSess = nSess - 2;
       
    % One or more sessions
    if isfield(aap.acq_details.input, 'referencedirectory_tmpl')
        refDir = strrep(aap.acq_details.input.referencedirectory_tmpl,'*',ID);
    end
    if ~isfield(aap.acq_details.input, 'selected_sessions') || ~any(aap.acq_details.input.selected_sessions),...
            aap.acq_details.input.selected_sessions = 1:nSess; 
    end
    
    if aap.directory_conventions.subject_directory_format == 3
        aap=aas_addsubject(aap,ID,VOL,'functional',aSess(aap.acq_details.input.selected_sessions));
    else
        aap=aas_addsubject(aap,VOL,'functional',aSess(aap.acq_details.input.selected_sessions));
    end
    
    for i = aap.acq_details.input.selected_sessions
		if ~aSess(i), continue; end
        
        % Obtain TR
        h = spm_dicom_headers(mri_finddcm(aap,VOL,aSess(i)));
        TR = h{1}.RepetitionTime/1000; % in seconds
        
        session = list_index(LIST{1},head.FMRI1,i+1);
        aap = aas_addsession(aap,session);
        if any(stagenumModel) && exist('refDir','var')
			clear names durations onsets tmod pmod;
            load(fullfile(refDir,['condition_vol_' ID '-' session '.mat']));
            % [MDV] initialise pmod and tmod if they don't exist
            if ~exist('pmod','var'), pmod(1:numel(names)) = struct('name',[],'param',[],'poly',[]);end
            if ~exist('tmod','var'), tmod = cell(1,numel(names)); end
            for iEV = 1:numel(names)
                % Event onsets has to be corrected according to the number of dummies
                % [MDV] format parametric from pmod and tmod if they are there
                if ~isempty(tmod{iEV}), % put any tmod in pmod
                    pmod(iEV).name = [{'time'} pmod(iEV).name];
                    pmod(iEV).param = [{onsets{iEV}-numdummies*TR} pmod(iEV).param];
                    pmod(iEV).poly = [tmod(iEV) pmod(iEV).poly];
                end
                l = numel(pmod(iEV).name);
                if l > 0,
                    clear parametric
                    parametric(1:l) = struct('name',[],'P',[],'h',[]);
                    for n = 1:l,
                        parametric(n).name = pmod(iEV).name{n};
                        parametric(n).P = pmod(iEV).param{n}';
                        parametric(n).h = pmod(iEV).poly{n};
                    end
                    aap=aas_addevent(aap,firstlevel,aap.acq_details.subjects(end).subjname,session,names{iEV},onsets{iEV}-numdummies*TR,durations{iEV},parametric);
                else
                    aap=aas_addevent(aap,firstlevel,aap.acq_details.subjects(end).subjname,session,names{iEV},onsets{iEV}-numdummies*TR,durations{iEV});
                end
            end
        end
    end
end
end

function ind = study_header(LIST)
LIST = importdata(LIST);
for n = 1:size(LIST,1)
    if list_search(LIST{n},'ID') % ID must be included
        i = 1;
        while ~isempty(list_index(LIST{n},i))
            ind.(list_index(LIST{n},i)) = i;
            i = i + 1;
        end
    end
end
end