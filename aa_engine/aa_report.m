% Automatic analysis - produces HTML summary of analysis
% This feature is usable, but still under development.
% Rhodri Cusack MRC CBU Cambridge 2004
% Tibor Auer MRC CBU Cambridge 2012-2013

function [aap]=aa_report_new(studyroot,stages)

if exist('studyroot','var')
    cd(studyroot);
else
    studyroot=pwd;
end;

% First, load AAP structure
load('aap_parameters');

if ~exist('stages','var')
    stages={aap.tasklist.main.module.name};
end;

if isfield(aap,'report'), aap = rmfield(aap,'report'); end
aap.report.numdependencies=0;

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
nSessions = numel(aap.acq_details.sessions);

for k=1:numel(stages)
    % domain
    xml = xml_read([stages{k} '.xml']);
    mfile_alias = stages{k};
    if ~exist(mfile_alias,'file'), mfile_alias = xml.tasklist.currenttask.ATTRIBUTE.mfile_alias; end
    domain = xml.tasklist.currenttask.ATTRIBUTE.domain;
    
    % Set inSession flag
    if strcmp(domain,'session') && ~inSession
        inSession = true;
    end
    if ~strcmp(domain,'session') && inSession
        inSession = false;
    end
    
    % Switch for stage
    all_stage = cell_index(stages, stages{k});
    istage = cell_index({aap.tasklist.main.module.name}, stages{k});
    istage = istage(all_stage==k);
    aapreport = aap.report;
    aap = aas_setcurrenttask(aap,istage);
    aap.report = aapreport;    
    
    % build dependency
    dep = aas_dependencytree_allfromtrunk(aap,domain);
	
	% run through 
    for d = 1:numel(dep)
        indices = dep{d}{2};
        isdone = exist(aas_doneflag_getpath_bydomain(aap,domain,indices,k),'file');
        try subj = dep{d}{2}(1); catch, subj = []; end % Subjects No
        try sess = dep{d}{2}(2); catch, sess = 1; end % Session/Occurrance No
        
        if sess == 1, aap = aas_report_add(aap,subj,['<h2>Stage: ' stages{k} '</h2>']); end
        
		% evaluate with handling sessions
        if inSession
            if sess == 1, aap = aas_report_add(aap,subj,'<table><tr>'); end % Open session
            aap = aas_report_add(aap,subj,'<td>');
            aap = aas_report_add(aap,subj,['<h3>Session: ' aap.acq_details.sessions(dep{d}{2}(2)).name '</h3>']);
        end
        if ~isdone
            aap = aas_report_add(aap,subj,'<h3>Not finished yet!</h3>');
        else
            [aap,resp]=aa_feval_withindices(mfile_alias,aap,'report',indices);
        end;
        if inSession
            aap = aas_report_add(aap,subj,'</td>');
            if sess == nSessions, aap = aas_report_add(aap,subj,'</tr></table>'); end % Close session
        end
        
    end
end

% Close files
aap = aas_report_add(aap,[],'EOF');
aap = aas_report_add(aap,0,'EOF');
fclose all;

% Show report
web(['file://' aap.report.html_main.fname]);
% Last, save AAP structure
save('aap_parameters_reported.mat', 'aap');
end