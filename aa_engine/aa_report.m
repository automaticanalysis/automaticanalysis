% Automatic analysis - produces HTML summary of analysis
%
% For reporting a remote pipeline, copy the aap_parameters.mat file to a
% local directory and call aa_report from that directory without any argument
%
% Rhodri Cusack MRC CBU Cambridge 2004
% Tibor Auer MRC CBU Cambridge 2012-2016

function [aap]=aa_report(studyroot,stages)

fprintf('Fetching report started...\n');

% First, load AAP structure
if exist('studyroot','var')
    if isstruct(studyroot)
        studyroot=fullfile(aas_getstudypath(studyroot),studyroot.directory_conventions.analysisid);
    end;
    cd(studyroot);
else
    studyroot=pwd;
end;
if ~exist('aap_parameters.mat','file')
    error('ERROR: aap structure not found');
else
    load('aap_parameters');
end

aa_init(aap);

if ~exist('stages','var')
    stages={aap.tasklist.main.module.name};
end;

if isfield(aap,'report'), aap = rmfield(aap,'report'); end
aap.report.numdependencies=0;
aap.report.dependency={};

aap.internal.total=0;
aap.internal.stagesnotdone=0;
stage_study_done = false;

% Provenance
aap.prov = aa_provenance(aap);

% Main HTLMs
if (isfield(aap.directory_conventions,'reportname'))
    aap.report.html_main.fname=fullfile(studyroot,aap.directory_conventions.reportname);
else
    aap.report.html_main.fname=fullfile(studyroot,'report.htm');
end;
aap.report.html_S00.fname = strrep(aap.report.html_main.fname,'.htm','_subjects.htm');
aap.report.subdir = fullfile(fileparts(aap.report.html_main.fname),'report_subjects');
aap.report.html_moco.fname = strrep(aap.report.html_main.fname,'.htm','_moco.htm');
aap.report.html_reg.fname = strrep(aap.report.html_main.fname,'.htm','_reg.htm');
aap.report.html_C00.fname = strrep(aap.report.html_main.fname,'.htm','_scon.htm');
aap.report.condir = fullfile(fileparts(aap.report.html_main.fname),'report_scon');

% Initialize HTMLs
aap = aas_report_add(aap,[],'HEAD=AA Report');
aap = aas_report_add(aap,0,'HEAD=Subjects');
if ~exist(aap.report.subdir,'dir'), mkdir(aap.report.subdir); end
aap = aas_report_add(aap,'moco','HEAD=Motion corection summary');
aap = aas_report_add(aap,'reg','HEAD=Registration summary');
aap = aas_report_add(aap,'C00','HEAD=First level contrasts');
if ~exist(aap.report.condir,'dir'), mkdir(aap.report.condir); end

aap.report.fbase = basename(aap.report.html_main.fname);

% Handle session
inSession = false;

for k=1:numel(stages)
    fprintf('Fetching report for %s...\n',stages{k});
    % domain
    if isfield(aap.tasklist.main.module(k),'aliasfor') && ~isempty(aap.tasklist.main.module(k).aliasfor)
        xml = xml_read([aap.tasklist.main.module(k).aliasfor '.xml']);
    else
        xml = xml_read([stages{k} '.xml']);
    end
    mfile_alias = stages{k};
    if ~exist(mfile_alias,'file'), mfile_alias = aap.tasklist.main.module(k).aliasfor; end
    if ~exist(mfile_alias,'file'), mfile_alias = xml.tasklist.currenttask.ATTRIBUTE.mfile_alias; end
    domain = xml.tasklist.currenttask.ATTRIBUTE.domain;
    
    % Switch for stage
    all_stage = cell_index(stages, stages{k});
    istage = cell_index({aap.tasklist.main.module.name}, stages{k});
    istage = istage(all_stage==k);
    
    % add provenance
    if aap.prov.isvalid, aap.prov.addModule(istage); end
    
    % Skip stages of unknown domain - most like aamod_importfilesasstream
    if ~strcmp(domain,'[unknown]')
        
        % build dependency
        dep = aas_dependencytree_allfromtrunk(aap,domain);
        
        % Set inSession flag
        if ~isempty(strfind(domain,'session')) && ~inSession
            inSession = true;
        end
        if isempty(strfind(domain,'session')) && inSession
            inSession = false;
        end
        
        % run through
        for d = 1:numel(dep)
            indices = dep{d}{2};
            isdone = exist(aas_doneflag_getpath_bydomain(aap,domain,indices,k),'file');
            if numel(dep{d}{2}), subj = dep{d}{2}(1); else subj = []; end % Subjects No
            if numel(dep{d}{2}) >= 2, sess = dep{d}{2}(2); end % Session/Occurrance No
            
            if inSession
                [junk, iSess] = aas_getN_bydomain(aap,domain,subj);
                firstSess = iSess(1);
                lastSess = iSess(end);
            end
            
            if ~inSession || (sess == firstSess), aap = aas_report_add(aap,subj,...
                    ['<h2>Stage: ' stages{k} aap.tasklist.main.module(k).extraparameters.aap.directory_conventions.analysisid_suffix '</h2>']); 
            end
            
            % evaluate with handling sessions
            if inSession
                if sess == firstSess, aap = aas_report_add(aap,subj,'<table><tr>'); end % Open session
                aap = aas_report_add(aap,subj,'<td valign="top">');
                aap = aas_report_add(aap,subj,['<h3>Session: ' aap.acq_details.([domain 's'])(dep{d}{2}(2)).name '</h3>']);
            end
            if ~isdone
                aap = aas_report_add(aap,subj,'<h3>Not finished yet!</h3>');
            else
                % set local aap
                if ~isempty(subj)
                    curr_aap = aas_setcurrenttask(aap,k,'subject',subj);
                else
                    curr_aap = aas_setcurrenttask(aap,k);
                end
                
                % set fields
                curr_aap.report = aap.report;
                curr_aap.prov = aap.prov;
         
                [curr_aap,junk]=aa_feval_withindices(mfile_alias,curr_aap,'report',indices);
                
                % save report
                aap.report = curr_aap.report;
                aap.prov = curr_aap.prov;
            end;
            if inSession
                aap = aas_report_add(aap,subj,'</td>');
                if sess == lastSess, aap = aas_report_add(aap,subj,'</tr></table>'); end % Close session
            end
        end
    end
end

% Close files
aap = aas_report_add(aap,[],'EOF');
aap = aas_report_add(aap,0,'EOF');
fclose all;

% Provenance
aap.prov.serialise(studyroot);

% Show report
web(['file://' aap.report.html_main.fname]);
% Last, save AAP structure
save(fullfile(studyroot,'aap_parameters_reported.mat'), 'aap');

aa_close;

end