% Automatic analysis - produces HTML summary of analysis
% This feature is usable, but still under development.
% Rhodri Cusack MRC CBU Cambridge 2004

function [aap]=aa_report(studyroot,stages)

if (~exist('studyroot','var'))
    studyroot=pwd;
end;


% First, load AAP structure
aaploadfn=fullfile(studyroot,'aap_parameters');
load(aaploadfn);

if (~exist('stages'))
    stages=aap.tasklist.stages;
end;

if (~iscell(stages))
    stages={stages};
end;

aap.report.html='';
aap.report.numdependencies=0;
aap.report.dependency={};

% make full report by calling individual components
% stage by stage display

aap.internal.total=0;
aap.internal.stagesnotdone=0;

aap.report.html=strcat(aap.report.html,'<table>');

for i=1:length(aap.acq_details.subjects)
    aap.report.html=strcat(aap.report.html,'<tr><td>');
    aap.report.html=strcat(aap.report.html,['<h1>Subject: ' aap.acq_details.subjects{i} '</h1>']);
    for k=1:length(stages)
        aap.report.html=strcat(aap.report.html,'<tr><td>');
        aap.report.html=strcat(aap.report.html,['<h2>Stage: ' stages{k}.name '</h2>']);

        doneflag=fullfile(aap.acq_details.root,['done_' stages{k}.name]);

        % find out whether this module needs to be executed once per study, subject or session
        [aap,domain]=feval(stages{k}.name,aap,'domain');

        aap.tasklist.currenttask.epiprefix=stages{k}.epiprefix;

        switch (domain)
            case 'study'
                doneflag=fullfile(aap.acq_details.root,['done_' stages{k}.name]);
                if (length(dir(doneflag)))
                    [aap,resp]=feval(stages{k}.name,aap,'report');
                end;
            case 'subject'
                doneflag=fullfile(aas_getsubjpath(aap,i),['done_' stages{k}.name]);
                if (length(dir(doneflag)))
                    [aap,resp]=feval(stages{k}.name,aap,'report',i);
                end;

            case 'session'
                aap.report.html=strcat(aap.report.html,'<table><tr>');
                for j=aap.acq_details.selected_sessions
                    aap.report.html=strcat(aap.report.html,'<td>');
                    aap.report.html=strcat(aap.report.html,['<h3>Session: ' aap.acq_details.sessions{j} '</h3>']);
                    doneflag=fullfile(aas_getsesspath(aap,i,j),['done_' stages{k}.name]);
                    if (length(dir(doneflag)))
                        [aap,resp]=feval(stages{k}.name,aap,'report',i,j);
                    end;
                    aap.report.html=strcat(aap.report.html,'</td>');
                end;
                aap.report.html=strcat(aap.report.html,'</tr></table>');
            otherwise
                aas_log(aap,1,sprintf('Unknown domain %s associated with stage %s',aap.tasklist.domain{k},stages{k}.name));
        end;
        aap.report.html=strcat(aap.report.html,'</td></tr>');
    end;
    aap.report.html=strcat(aap.report.html,'</td></tr>');

end;
aap.report.html=strcat(aap.report.html,'</table>');

htm=aap.report.html;
if (isfield(aap.directory_conventions,'reportname'))
    reportname=aap.directory_conventions.reportname;
else
    reportname='fullreport.htm';
end;

fid=fopen(fullfile(aap.acq_details.root,reportname),'w');
fprintf(fid,'%s',htm);
fclose(fid);
myurl=['file://' fullfile(aap.acq_details.root,reportname)];
web(myurl);

return;

