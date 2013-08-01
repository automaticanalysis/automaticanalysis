% Automatic analysis - produces HTML summary of analysis
% This feature is usable, but still under development.
% Rhodri Cusack MRC CBU Cambridge 2004
% Tibor Auer MRC CBU Cambridge 2012-2013

function [aap]=aa_report(studyroot,stages)

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
aap.report.dependency={};

aap.internal.total=0;
aap.internal.stagesnotdone=0;
stage_study_done = false;

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

aap = aas_report_add(aap,'HEAD=AA Report');
aap = aas_report_add(aap,0,'HEAD=Subjects');
if ~exist(aap.report.subdir,'dir'), mkdir(aap.report.subdir); end
aap = aas_report_add(aap,'moco','HEAD=Motion corection summary');
aap = aas_report_add(aap,'reg','HEAD=Registration summary');
aap = aas_report_add(aap,'C00','HEAD=First level contrasts');
if ~exist(aap.report.condir,'dir'), mkdir(aap.report.condir); end

aap.report.fbase = basename(aap.report.html_main.fname);
for i=1:numel(aap.acq_details.subjects)
    fprintf('Reporting subject: %s\n',basename(aas_getsubjpath(aap,i)));
    aap.report.(sprintf('html_S%02d',i)).fname = fullfile(aap.report.subdir,[aap.report.fbase sprintf('_S%02d.htm',i)]);
    aap = aas_report_add(aap,0,...
        sprintf('<a href="%s" target=_top>%s</a><br>',...
        aap.report.(sprintf('html_S%02d',i)).fname,...
        ['Subject: ' basename(aas_getsubjpath(aap,i))]));
    
    % Single-subject reports
    aap = aas_report_add(aap,i,['HEAD=Subject: ' basename(aas_getsubjpath(aap,i))]);
    for k=1:numel(stages)
        
        % find out whether this module needs to be executed once per study, subject or session
        xml = xml_read([stages{k} '.xml']);
        mfile = stages{k};
        if ~exist(mfile,'file'), mfile = xml.tasklist.currenttask.ATTRIBUTE.mfile_alias; end
        domain = xml.tasklist.currenttask.ATTRIBUTE.domain;
        
        % handle repetition
        all_stage = cell_index(stages, stages{k});
        istage = cell_index({aap.tasklist.main.module.name}, stages{k});
        istage = istage(all_stage==k);
        
        switch (domain)
            case 'study'
                if ~stage_study_done
                    aap = aas_report_add(aap,['<h2>Stage: ' stages{k} '</h2>']);
                    doneflag=fullfile(aas_getstudypath(aap,k),['done_' aas_getstagetag(aap,k)]);
                    if ~exist(doneflag,'file')
                        aap = aas_report_add(aap,i,'<h3>Not finished yet!</h3>');
                    else
%                         curr_aap = load(fullfile(fileparts(aas_getsubjpath(aap,i,k)),['aap_parameters_' aas_getstagetag(aap,k) '.mat'])); curr_aap = curr_aap.aap;
                        aapreport = aap.report;
                        aap = aas_setcurrenttask(aap,istage);
                        aap.report = aapreport;
                        aap = feval(mfile,aap,'report');
                    end;
                    stage_study_done = true;
                end
            case 'subject'
                aap = aas_report_add(aap,i,['<h2>Stage: ' stages{k} '</h2>']);                
                doneflag=fullfile(aas_getsubjpath(aap,i,k),['done_' aas_getstagetag(aap,k)]);
                if ~exist(doneflag,'file')
                    aap = aas_report_add(aap,i,['<h3>Not finished yet!</h3>']);
                else
%                     curr_aap = load(fullfile(aas_getsubjpath(aap,i,k),['aap_parameters_' aas_getstagetag(aap,k) '.mat'])); curr_aap = curr_aap.aap;
%                     curr_aap.acq_details.subjects = aap.acq_details.subjects;
                    aapreport = aap.report;
                    aap = aas_setcurrenttask(aap,istage);
                    aap.report = aapreport;
                    aap = feval(mfile,aap,'report',i);
                end;
            case 'session'
                aap = aas_report_add(aap,i,['<h2>Stage: ' stages{k} '</h2>']);                                
                aap = aas_report_add(aap,i,'<table><tr>');
                for j=1:numel(aap.acq_details.sessions)
                    aap = aas_report_add(aap,i,'<td>');
                    aap = aas_report_add(aap,i,['<h3>Session: ' aap.acq_details.sessions(j).name '</h3>']);
                    doneflag=fullfile(aas_getsesspath(aap,i,j,k),['done_' aas_getstagetag(aap,k)]);
                    if ~exist(doneflag,'file')
                        aap = aas_report_add(aap,i,['<h3>Not finished yet!</h3>']);
                    else
%                         curr_aap = load(fullfile(aas_getsesspath(aap,i,j,k),['aap_parameters_' aas_getstagetag(aap,k) '.mat'])); curr_aap = curr_aap.aap;
%                         curr_aap.acq_details.subjects = aap.acq_details.subjects;                        
                        aapreport = aap.report;
                        aap = aas_setcurrenttask(aap,istage);
                        aap.report = aapreport;
                        aap = feval(mfile,aap,'report',i,j);
                    end;
                    aap = aas_report_add(aap,i,'</td>');
                end;
                aap = aas_report_add(aap,i,'</tr></table>');
            otherwise
                aas_log(aap,1,sprintf('Unknown domain %s associated with stage %s',aap.tasklist.domain{k},stages{k}));
        end;
    end;
    aap = aas_report_add(aap,i,'EOF');
end;

aap = aas_report_add(aap,'EOF');

aap = aas_report_add(aap,0,'EOF');

fclose all;
web(['file://' aap.report.html_main.fname]);
% Last, save AAP structure
save('aap_parameters_reported.mat', 'aap');

end