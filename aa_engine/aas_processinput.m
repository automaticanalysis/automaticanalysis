% Allow specify Subjects, Sessions and Events based on text file.
% Required options in aap.acq_details.input (see also aap_parameters_defaults.xml):
%	- list: text file:
%		Required format:
%		- a CSV file: cells are separated with semicolon
%		- subcells are separated with "_"
%		- a header in the first line
%	 	Required columns (may contain more): 
%		- "ID": numbers for identifying subjects
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

% Reads in header
head = study_header(aap.acq_details.input.list);
LIST = importdata(aap.acq_details.input.list);
for v = 2:size(LIST,1)
    ID = str2double(list_index(LIST{v},head.ID));
    VOL = list_index(LIST{v},head.FMRI1);
    nSess = 2; aSess = [];
    while ~isempty(list_index(LIST{1},head.FMRI1,nSess))
        aSess = horzcat(aSess,str2double(list_index(LIST{v},head.FMRI1,nSess)));
        nSess = nSess + 1;
    end
    nSess = nSess - 2;
    strSubj = mri_findvol(str2double(VOL));
    
    % Obtain TR from the first session
    h = dicominfo(mri_finddcm(str2double(VOL),aSess(1)));
    TR = h.RepetitionTime/1000; % in seconds
    
    % One or more sessions
    if isfield(aap.acq_details.input, 'referencedirectory_tmpl')
        refDir = strrep(aap.acq_details.input.referencedirectory_tmpl,'*',num2str(ID));
    end
    if ~isfield(aap.acq_details.input, 'selected_sessions') || ~any(aap.acq_details.input.selected_sessions),...
            aap.acq_details.input.selected_sessions = 1:nSess; 
    end
    aap=aas_addsubject(aap,strSubj,aSess(aap.acq_details.input.selected_sessions));
    for i = aap.acq_details.input.selected_sessions
        session = list_index(LIST{1},head.FMRI1,i+1);
        aap = aas_addsession(aap,session);
        if exist('refDir','var')
            load(fullfile(refDir,['condition_vol_' num2str(ID) '-' session '.mat']));
            for iEV = 1:numel(names)
                % Event onsets has to be corrected accoring to the number of dummies
                aap=aas_addevent(aap,'aamod_firstlevel_model',strSubj,session,names{iEV},onsets{iEV}-aap.acq_details.numdummies*TR,durations{iEV});
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